# Variant calling

import collections
import numpy as np

import pavlib

DEFAULT_MIN_ANCHOR_SCORE = 1000


def call_from_align(caller_resources, min_anchor_score=DEFAULT_MIN_ANCHOR_SCORE, verbose=False):
    """
    Create a list of variant calls from alignment table.

    :param caller_resources: Caller resources.
    :param min_anchor_score: Minimum allowed score for an alignment segment to anchor a variant call.
    :param verbose: Print progress.

    :return: A list of variant call objects.
    """

    variant_call_list = list()

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
            if start_index not in qryref_index_set:
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
        # Resolve complex loci
        #

        sv_dict = dict()  # Key: interval range (tuple), value=SV object

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

            # Complete candidate variant call
            sv_dict[(chain_node.start_index, chain_node.end_index)] = variant_call


        #
        # Choose optimal SVs
        #

        # Initialize Bellman-Ford
        top_score = np.full(df_align.shape[0], -np.inf)  # Score, top-sorted graph
        top_tb = np.full(df_align.shape[0], -2)      # Traceback (points to parent node with the best score), top-sorted graph

        for source_node in chain_container.source_node_list:
            top_score[source_node.start_index] = df_align.loc[source_node.start_index]['SCORE']
            top_tb[source_node.start_index] = -1  # Points to root node

        # if np.isneginf(top_score[0]):
        #     #raise RuntimeError('Root node does not point to node 0')
        #     print('Root node does not point to node 0')
        #     continue

        # Create a graph by nodes (anchor graph nodes are scored edges)
        node_link = collections.defaultdict(list)

        for start_index, end_index in chain_container.chain_dict.keys():
            node_link[start_index].append(end_index)

        # Update score by Bellman-Ford
        for start_index in range(df_align.shape[0]):
            base_score = top_score[start_index]

            if np.isneginf(base_score):  # Unreachable
                continue

            for end_index in node_link[start_index]:
                score = base_score + sv_dict[start_index, end_index].score_variant  # Base + variant score
                score += caller_resources.score_model.gap(df_align.loc[end_index, 'QRY_POS'] - df_align.loc[start_index, 'QRY_END'])  # Unaligned query bases between records

                if score > top_score[end_index]:
                    top_score[end_index] = score
                    top_tb[end_index] = start_index

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

        variant_call_list.extend([sv_dict[node_interval] for node_interval in optimal_interval_list])

    return variant_call_list
