import os
import sys

def run( target_list , t_type):
    #'''
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # modify the deftype.h file
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!
    f = open( 'typedef_temp.h' , 'r')
    content = f.readlines() 
    print target_list 
    for lines in range(len(content)):
        lines_split = content[lines].split(' ')
        for target in target_list:
            if len(lines_split) >= 3 and lines_split[2] == target :
                print lines_split
                lines_split[1] = t_type 
        for idx in range(len(lines_split)-1):
            lines_split[idx] += ' '
        content[lines] = "".join(lines_split)
    
    content = ''.join(content)
    f.close()
    g = open( 'typedef.h' , 'w')
    g.write(content)
    g.close()
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #'''
    os.system('rm ./lulesh')
    os.system('rm temp_output.csv')
    os.system('make')
    os.system('./lulesh -s 50 -f 0 1 ')
