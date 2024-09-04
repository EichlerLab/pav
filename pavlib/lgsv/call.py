# Variant calling

import collections
import numpy as np

import pavlib
import svpoplib

DEFAULT_MIN_ANCHOR_SCORE = 1000


def call_from_align(caller_resources, min_anchor_score=DEFAULT_MIN_ANCHOR_SCORE, verbose=False, dot_filename=None):
    """
    Create a list of variant calls from alignment table.

    :param caller_resources: Caller resources.
    :param min_anchor_score: Minimum allowed score for an alignment segment to anchor a variant call.
    :param verbose: Print progress.
    :param dot_filename: Name of dot graph file to write if not None. File is gzipped if ends with the filename ends
        with '.gz'.

    :return: A list of variant call objects.
    """

    variant_call_list = list()
    variant_id_set = set()

    for query_id in caller_resources.df_align_qry['QRY_ID'].unique():
        if verbose:
            print(f'Query: {query_id}')

        df_align = caller_resources.df_align_qry.loc[
            caller_resources.df_align_qry['QRY_ID'] == query_id
        ].reset_index(drop=True)

        chain_container = pavlib.lgsv.chain.AnchorChainContainer()

        start_index = 0
        last_index = df_align.shape[0]

        qryref_index_set = set(caller_resources.df_align_qryref.index)

        # Traverse interval setting each position to the left-most anchor candidate.
        while start_index < last_index:

            # Skip if anchor did not pass TIG & REF trimming
            if df_align.loc[start_index]['INDEX'] not in qryref_index_set:
                start_index += 1
                continue

            start_row = df_align.loc[start_index]
            end_index = start_index + 1

            # Traverse each interval after the start for right-most anchor candidates. Limit search by contig distance
            while end_index < last_index and pavlib.lgsv.chain.can_reach_anchor(start_row, df_align.loc[end_index], caller_resources.score_model):

                if df_align.loc[end_index]['INDEX'] in qryref_index_set and pavlib.lgsv.chain.can_anchor(
                        start_row, df_align.loc[end_index], caller_resources.score_model, min_anchor_score
                ):
                    chain_container.add_anchor(start_index, end_index)

                end_index += 1

            start_index += 1


        #
        # Resolve variant calls
        #

        sv_dict = dict()  # Key: interval range (tuple), value=SV object
        min_sv_score = 0

        for chain_node in chain_container.chain_dict.values():

            if (chain_node.start_index, chain_node.end_index) in sv_dict:
                raise RuntimeError(f'Duplicated chain entries for index: {chain_node}')

            #
            # SV call fundamentals
            #
            interval = pavlib.lgsv.interval.AnchoredInterval(chain_node, df_align, caller_resources.score_model)

            # Try INS
            variant_call = pavlib.lgsv.variant.InsertionVariant(interval, caller_resources)

            # Try DEL
            variant_call_next = pavlib.lgsv.variant.DeletionVariant(interval, caller_resources)

            if variant_call_next.score_variant > variant_call.score_variant:
                variant_call = variant_call_next

            # Try Inversion
            variant_call_next = pavlib.lgsv.variant.InversionVariant(interval, caller_resources)

            if variant_call_next.score_variant > variant_call.score_variant:
                variant_call = variant_call_next

            # Try tandem duplication
            variant_call_next = pavlib.lgsv.variant.TandemDuplicationVariant(interval, caller_resources)

            if variant_call_next.score_variant > variant_call.score_variant:
                variant_call = variant_call_next

            # Try complex
            variant_call_next = pavlib.lgsv.variant.ComplexVariant(interval, caller_resources)

            if variant_call_next.score_variant > variant_call.score_variant:
                variant_call = variant_call_next

            # Set to Null variant if anchors cannot support the variant call
            if not variant_call.is_null() and variant_call.score_variant + variant_call.anchor_score_min < 0:
                variant_call = pavlib.lgsv.variant.NullVariant()

            if not variant_call.is_null():
                min_sv_score = min(min_sv_score, variant_call.score_variant)

            # Complete candidate variant call
            sv_dict[(chain_node.start_index, chain_node.end_index)] = variant_call


        #
        # Choose optimal SVs
        #

        # Initialize Bellman-Ford
        top_score = np.full(df_align.shape[0], -np.inf)   # Score, top-sorted graph
        last_aligned = np.full(df_align.shape[0], False)  # True if the last node was aligned
        top_tb = np.full(df_align.shape[0], -2)           # Traceback (points to parent node with the best score), top-sorted graph

        for source_node in chain_container.source_node_list:
            top_score[source_node.start_index] = 0
            top_tb[source_node.start_index] = -1  # Points to root node

        # if np.isneginf(top_score[0]):
        #     #raise RuntimeError('Root node does not point to node 0')
        #     print('Root node does not point to node 0')
        #     continue

        # Create a graph by nodes (anchor graph nodes are scored edges)
        node_link = collections.defaultdict(set)

        for start_index, end_index in chain_container.chain_dict.keys():
            node_link[start_index].add(end_index)

        for start_index in range(df_align.shape[0] - 1):
            node_link[start_index].add(start_index + 1)

        # Update score by Bellman-Ford
        for start_index in range(df_align.shape[0]):
            base_score = top_score[start_index]

            if np.isneginf(base_score):  # Unreachable
                continue

            for end_index in sorted(node_link[start_index]):
                sv_score = sv_dict[start_index, end_index].score_variant \
                    if (start_index, end_index) in sv_dict else -np.inf

                if not np.isneginf(sv_score):
                    score = base_score + \
                            sv_score + \
                            df_align.loc[end_index]['SCORE'] + \
                            df_align.loc[start_index]['SCORE'] if not last_aligned[start_index] else 0  # Add anchor score for left anchor if not in the chain
                    #caller_resources.score_model.gap(df_align.loc[end_index, 'QRY_POS'] - df_align.loc[start_index, 'QRY_END'])  # Unaligned query bases between records
                    #sv_score = min_sv_score * (end_index - start_index)
                else:
                    score = base_score  # Skip over null variants

                #score = base_score + sv_dict[start_index, end_index].score_variant  # Base + variant score
                #score = base_score + sv_score + df_align.loc[end_index]['SCORE']  # Base + variant score + right anchor
                #score += caller_resources.score_model.gap(df_align.loc[end_index, 'QRY_POS'] - df_align.loc[start_index, 'QRY_END'])  # Unaligned query bases between records

                if score > top_score[end_index]:
                    top_score[end_index] = score
                    top_tb[end_index] = start_index
                    last_aligned[end_index] = not np.isneginf(sv_score)

        # Trace back through scored paths
        sink_node_list = sorted({
            chain_node.end_index
                for chain_node in chain_container.chain_dict.values()
                    if len(chain_node.next_node_list) == 0
        })

        if len(sink_node_list) == 0:
            continue

        max_score_index = np.argmax([top_score[i] for i in sink_node_list])
        max_sink_node = sink_node_list[max_score_index]

        optimal_interval_list = list()

        last_node = max_sink_node

        while True:
            first_node = top_tb[last_node]

            if first_node < 0:
                break

            optimal_interval_list.append((first_node, last_node))

            last_node = first_node

        # Add variants and version IDs
        new_var_list = [
            sv_dict[node_interval] for node_interval in optimal_interval_list \
                if node_interval in sv_dict and not sv_dict[node_interval].is_null()
        ]

        for variant_call in new_var_list:
            variant_call.variant_id = svpoplib.variant.version_id_name(variant_call.variant_id, variant_id_set)
            variant_id_set.add(variant_call.variant_id)

        variant_call_list.extend(new_var_list)

        # Write dot file
        if not dot_filename is None:
            with pavlib.io.PlainOrGzFile(dot_filename, 'wt') as out_file:
                pavlib.lgsv.util.dot_graph_writer(
                    out_file, df_align, chain_container, sv_dict, graph_name=f'"{query_id}"', forcelabels=True
                )

    # Return list
    return variant_call_list
