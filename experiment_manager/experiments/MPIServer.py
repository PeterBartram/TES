# TES is an open source integration package for modelling exoplanet evolution.
# Copyright (C) <2021>  <Peter Bartram, Alexander Wittig>

# TES is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>. 

import numpy as np
from mpi4py import MPI as mpi
import subprocess
import copy 

def running_single_core():
    comm = mpi.COMM_WORLD
    size = comm.Get_size()

    if size == 1:
        return True
    else:
        return False
    
def distribute_seed(seed=0):
    data = None
    comm = mpi.COMM_WORLD
    size = comm.Get_size()
    
    if size > 1:
        rank = comm.Get_rank()
        
        if rank == 0:
            data = {}
            data['seed'] = seed
            
        data = comm.bcast(data)
        return data['seed']
    else:
        return seed


def run_experiments(experiments):       
    exp_count = len(experiments)
    comm = mpi.COMM_WORLD
    comm.Barrier()
    size = comm.Get_size()
    rank = comm.Get_rank()
    
    if rank == 0:       
          sendCount = 0
          recCount = 0
      
          initial_burst_size = np.minimum(size, exp_count+1)

          # Distribute data across MPI nodes. 
          for i in range(1, initial_burst_size):
              message = {}
              message['experiment index'] = sendCount
              message['terminate'] = False
              request = comm.isend(message, dest=i, tag=1)    
              request.Wait()
              sendCount += 1
        
          # Terminate any threads that are not in use. 
          for i in range(initial_burst_size, size):
              message = {}
              message['terminate'] = True
              request = comm.isend(message, dest=i)          
              request.Wait()
              
          # Data sever - send out experiment numbers to be run. 
          while(recCount < exp_count):
              message = comm.recv()
              recCount += 1
              source = int(message['source'])
      
              # Index of experiment performed
              index = np.int(message['experiment index'])
             
              # Copy the results back to the server experiment database.
              experiments[index].integration = message['integration']
              
              if (sendCount >= exp_count):
                  print("Sending termination to:", source)
                  # Send message to terminate this thread.
                  message['terminate'] = True
                  message['count'] = -1
                  request = comm.isend(message, source)          
                  request.Wait()
              else:
                  # Send message to start another integration.
                  message['terminate'] = False
                  message['experiment index'] = sendCount
                  sendCount += 1
                  request = comm.isend(message, source)                
                  request.Wait()
      
          print("Finished all integrations.", flush=True)
          return experiments
    
    else:
        # Data client - run the experiments and send the data back to the server.
        while(1):
            message = comm.recv()
            
            if message['terminate'] == True:
                print('Exiting on rank', rank, flush=True)
                # Finished with server MPI threads so kill them off.
                exit()
            else:
                print("Executing integration number", message['experiment index'])
        
                experiments[message['experiment index']].run()
                experiments[message['experiment index']].save()

                message['integration'] = copy.deepcopy(experiments[message['experiment index']].integration)
                message['source'] = rank
                request = comm.isend(message, dest=0, tag=1)
                print("Sent reply to thread 0.")        
                request.Wait()
                
        
        
