import multiprocessing
from typing import List, Callable
import numpy as np


class multiprocessed_dynamic_programming:
    #Segment scoring functions must return negative values since the class tries to maximize score
    scoring_functions: List[Callable]
    #Accumulating function dictates how scores of multiple segments are combined
    accumulating_function: Callable
    #Stored parameter list containing all necessary information
    parameters: dict

    #Helper function for multiprocessing
    def consumer(self, task_queue, finish_queue, task_object, task_name, results):
        #processes element in queue
        while(True):
            idx = task_queue.get()
            try:
                getattr(task_object, task_name)(idx, results)
                finish_queue.put((idx, False))
            except Exception as e:
                finish_queue.put((None, e))  # Indicates an error occurred
                break  # Exit the loop if an error occurs
    
    def run_singleprocess_calculation(self, length, types, protection_length, task_object, task_name, verbose = None, coarseness = 1):
        #results = [multiprocessing.Array(array_type, length) for array_type in types]
        results = [np.zeros(length, dtype = array_type) for array_type in types]
        for results_array in results:
            for i in range(len(results_array)):
                results_array[i] = -1
        for i in range(0, length, coarseness):
            getattr(task_object, task_name)(i, results)
            if verbose is not None:
                verbose(i, length)
        getattr(task_object, task_name)(length - 1, results)
        return results

    def run_multiprocess_calculation(self, num_workers, length, types, protection_length, task_object, task_name, verbose, coarseness):
        assert (protection_length > 0)
        self.protection_length = (protection_length // coarseness) * coarseness
        if self.protection_length < protection_length:
            self.protection_length = self.protection_length + coarseness
        if(num_workers == 1):
            return self.run_singleprocess_calculation(self, length, types, protection_length, task_object, task_name, verbose, coarseness)
        #Initialize Multiprocessing Queues
        task_queue = multiprocessing.Queue()
        finish_queue = multiprocessing.Queue()

        #Created Shared Array
        results = [multiprocessing.Array(array_type, length, lock = False) for array_type in types]
        for results_array in results:
            for i in range(len(results_array)):
                results_array[i] = -1

        #Create Worker Processes
        workers = [multiprocessing.Process(target = self.consumer, args = (self, task_queue, finish_queue, task_object, task_name, results)) for _ in range(num_workers)]

        #Start Worker Processes
        for worker in workers:
            worker.start()

        #Initialize list to keep track of completed tasks
        processed = [False]*length

        #Task Creation
        for i in range(0, self.protection_length, coarseness):
            task_queue.put(i)

        min_unfinished_index = 0
        #Minimum index such that all indices less than min_unfinished_index are True

        max_index = (length // coarseness) * coarseness

        while min_unfinished_index < length:
            # Collect New Finished Task
            idx, error = finish_queue.get()
            if error:
                print("Error in Worker:")
                # Terminate all workers and join
                for worker in workers:
                    worker.terminate()
                    worker.join()
                raise Exception(error)
            processed[idx] = True

            # Update min_unfinished_index and create new tasks
            while min_unfinished_index < length and processed[min_unfinished_index]:
                if min_unfinished_index + self.protection_length < length: # Create new task
                    task_queue.put(min_unfinished_index + self.protection_length)
                if (verbose is not None):
                    verbose(min_unfinished_index, length)
                min_unfinished_index = min_unfinished_index + coarseness
            if (min_unfinished_index > length and not processed[length - 1]):
                min_unfinished_index = length - 1
                task_queue.put(min_unfinished_index)
        #Get rid of worker
        for worker in workers:
            worker.terminate()
            worker.join()
        return [list(result) for result in results]