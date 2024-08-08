for i in *.avi; do ffmpeg -i "$i" -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" "${i%.*}.mp4"; done
# for i in *.mp4; do ffmpeg -i "$i" -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" "${i%.*}_2.mp4"; done
read -p $'\e[33mRemove the .avi files here? [-1]:\e[0m' remove_avi_files
remove_avi_files=${remove_avi_files:-1}
if [[ remove_avi_files == 1 ]]; then
  rm *.avi
