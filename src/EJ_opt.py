from scipy.optimize import differential_evolution
import EJ_Run
import gen_targets

t_type   = 'float'
baseline_e = 512427 
#target_candid = gen_targets.CalcKine_kerenl
#target_candid = gen_targets.CalcVolumeFor_kernel
#target_candid = gen_targets.ApplyMPAndUV_kernel
target_candid = gen_targets.all_targets
Opt_results = open('Opt_results.csv' , 'a')
Opt_results.write( '=========================================================\n')

def Lulush_out( target_idx ):
    target_list = get_target( target_idx)
    EJ_Run.run( target_list , t_type );
    is_Error = False 
    try:
        f = open('temp_output.csv' , 'r' ) 
        result = f.readlines()
        f.close()
        print result 
        result_sp = result[0][0:-1].split(',')
        runTime = float(result_sp[0])
        Energy =  float(result_sp[1])
    except:
        print '***** Compile Error*******'
        is_Error = True
        runTime = float('inf')
        Energy =  float('inf')
    
    
    for targets in target_list:
        Opt_results.write(targets)
        Opt_results.write(' ')
    Opt_results.write( ',')
    Opt_results.write( str(is_Error) )
    Opt_results.write( ',')
    Opt_results.write( str(runTime) )
    Opt_results.write( ',')
    Opt_results.write( str(Energy) )
    Opt_results.write( '\n')
    
    if abs(Energy - baseline_e )/baseline_e:
        return runTime
    else:
        return float('inf')
   
#def get_target( target_candid , idx ):
def get_target( idx ):
    target_list = []
    for isChoose in range(len(idx)):
        if idx[isChoose] > 0.5:
            target_list.append(target_candid[isChoose]) 
    #print target_list
    #return len(idx) - len(target_list)
    return target_list 

# Target_candid
tar_idx = []
for x in range(len(target_candid)):
    tar_idx.append((0,1.0))
#print get_target( tar_can  , tar_idx )
final_result = differential_evolution( Lulush_out , tar_idx )
print final_result
Opt_results.close()
#print Lulush_out( ['Real_t_x'] ) 
