for i in *.pdf; do
   convert           \
   -verbose       \
   -density 150   \
   -trim          \
    $i      \
   -quality 100   \
   -sharpen 0x1.0 \
    $i.png
done
