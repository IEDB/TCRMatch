#!/usr/bin/bash

if make && cd scripts && bash update.sh && cd ../; then
    echo "Install successful"
else
    echo "An error occurred during installation"
fi
