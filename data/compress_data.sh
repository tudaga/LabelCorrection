#!/usr/bin/env bash

echo "compressing data" 
for file in $(find . -name '*.csv' -or -name '*.h5ad'); do
	tarfile="${file%.*}.tar.lzma"
	if [[ -f "$tarfile" ]]; then
		echo "found ${tarfile}, skipping"
	else
		tar -cv --lzma -f $tarfile $file
		rm $file
		echo "compressed ${file}"

	fi
done

