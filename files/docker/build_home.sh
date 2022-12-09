# Build dot files in root and construct a default root for the run script.

# Set root
cat << EOF >> /etc/bash.bashrc

# System-wide aliases
alias ls='ls -aFh --color=auto'
alias ll='ls -aFlh --color=auto'
EOF

# Setup default
mkdir -p /home/default -m 777
