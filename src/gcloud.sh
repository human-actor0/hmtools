# here some useful gcloud tools are collected

gcloud.upload(){
usage="$FUNCNAME <input> <google bucket> [<options>]"
if [ $# -lt 2 ];then echo "$usage"; return; fi
	gsutil -o GSUtil:parallel_composite_upload_threshold=150M cp $@
}

