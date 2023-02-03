#!/usr/bin/env bash

for file in $(find . -name '*.tar.lzma'); do
	tar –xvzf $file
done

