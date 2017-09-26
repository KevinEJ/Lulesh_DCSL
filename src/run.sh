BASIC="./lulesh -s 50 -f"

rm output.csv

#for BLOCK_SIZE in 1 2 4 8 16 32 64 128 256 512 1024
for BLOCK_SIZE in 128 256 512 1024
do
    #for FUNC_NUM in 1 2 3 4 5 6 7 8
    for FUNC_NUM in 2 6 7
    do
        $BASIC $FUNC_NUM $BLOCK_SIZE
    done
done

