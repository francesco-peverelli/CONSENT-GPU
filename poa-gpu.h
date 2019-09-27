#ifndef POA_GPU_H
#define POA_GPU_H

#include <cuda_runtime.h>
#include <chrono>
#include <unistd.h>
#include <thrust/scan.h>
#include <thrust/device_vector.h>
#include <thrust/device_ptr.h>
#include <numeric>
#include <stdexcept>
#include <thread>
#include "cuda-poa.cuh"
#include "tasks_manager.h"

using namespace std;
using namespace chrono;

#define NOW high_resolution_clock::now() 

#define cudaErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char* file, int line, bool abort = true) {

	if (code != cudaSuccess) {
		fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
		if (abort) exit(code);
	}
}

namespace poa_gpu_utils{

template<int MAX_SEQ>
vector<string> form_window(char* sequences, const int nseq, const int seq_size) {
	
	char* seq_ptr = sequences;
	char tmp_seq[MAX_SEQ+1]; 
	
	vector<string> result(nseq);
	for(int j = 0; j < nseq; j++){
		
		for(int k = 0; k < MAX_SEQ; k++){
			tmp_seq[k] = seq_ptr[k];
		}
		//memcpy(tmp_seq, seq_ptr, MAX_SEQ_32);
		tmp_seq[seq_size] = '\0';
		result[j] = tmp_seq;
		seq_ptr += seq_size;
	}
	return result;
}

template<int BDIM>
void init_kernel_block_parameters(vector<Task<vector<string>>> &window_batch, char** sequences, vector<int> &nseq_offsets, vector<int> &seq_offsets, int* tot_nseq, int first_el) {
	
	int batch_size = window_batch.size();						
	int n = (first_el + BDIM < batch_size) ? BDIM : batch_size - first_el;
	for(int window_idx = first_el; window_idx < first_el + BDIM && window_idx < batch_size; window_idx++) {
		
		nseq_offsets[window_idx - first_el] = window_batch[window_idx].task_data.size();
	}
	partial_sum(nseq_offsets.begin(),nseq_offsets.end(),nseq_offsets.begin());
	*tot_nseq = nseq_offsets[n-1];
	seq_offsets = vector<int>(*tot_nseq);
	
	int sequence_idx = 0;
	for(int window_idx = first_el; window_idx < first_el + BDIM && window_idx < window_batch.size(); window_idx++) {
		vector<string> &window = window_batch[window_idx].task_data;			
		int wsize = window.size();
		for(int i = 0; i < wsize; i++) {
			seq_offsets[sequence_idx] = window[i].size();
			sequence_idx++;
		}
	}
	partial_sum(seq_offsets.begin(), seq_offsets.end(), seq_offsets.begin());
	int tot_size = seq_offsets[sequence_idx-1];
	*sequences = (char*)malloc(tot_size);
	sequence_idx = 0;

	for(int window_idx = first_el; window_idx < first_el + BDIM && window_idx < window_batch.size(); window_idx++) {
			
		vector<string>& window = window_batch[window_idx].task_data;	
		for(auto seq : window) {
			int offset;
			if(sequence_idx == 0){
				offset = 0;
			}else{
				offset = seq_offsets[sequence_idx-1];
			}
			char* seq_ptr = (*sequences) + offset;
			memcpy(seq_ptr, seq.c_str(), seq.size());
			sequence_idx++;
		}
	}	
}

template<int SL, int MAXL, int WL, int BDIM>                                                                                                      inline void gpu_POA_alloc(TaskRefs &T){

	//cout << "Allocating memory for " << SL << " " << WL << " " << MAXL << " " << BDIM << "\n";

	T.result = (char*)malloc(WL * MAXL * BDIM);
	T.res_size = (int*)malloc(BDIM * sizeof(int));

	cudaErrchk(cudaMalloc(&T.old_len_global_d, (unsigned long long)BDIM * sizeof(int)));
	cudaErrchk(cudaMalloc(&T.dyn_len_global_d, (unsigned long long)BDIM * sizeof(int)));

	cudaErrchk(cudaMalloc(&T.sequences_d, (unsigned long long)SL * WL * BDIM)); 
	cudaErrchk(cudaMalloc(&T.seq_offsets_d, (unsigned long long)BDIM * WL * sizeof(int))); 
	cudaErrchk(cudaMalloc(&T.nseq_offsets_d, (unsigned long long)BDIM * sizeof(int))); 
	cudaErrchk(cudaMalloc(&T.result_d, (unsigned long long)MAXL * WL * BDIM)); 
	cudaErrchk(cudaMalloc(&T.res_size_d, (unsigned long long)BDIM * sizeof(int)));

	cudaErrchk(cudaMalloc(&T.lpo_edge_offsets_d, (unsigned long long)WL * BDIM * sizeof(int)));
	cudaErrchk(cudaMalloc(&T.lpo_letters_d, (unsigned long long)WL * SL * BDIM));
	cudaErrchk(cudaMalloc(&T.lpo_edges_d, (unsigned long long)WL * SL * EDGE_F * BDIM * sizeof(Edge)));
	cudaErrchk(cudaMalloc(&T.edge_bounds_d, (unsigned long long)WL * (SL+1) * BDIM * sizeof(int)));	
	cudaErrchk(cudaMalloc(&T.end_nodes_d, (unsigned long long)WL * SL * BDIM));
	cudaErrchk(cudaMalloc(&T.sequence_ids_d, (unsigned long long)WL * WL * SL * BDIM));
	
	cudaErrchk(cudaMalloc(&T.new_letters_global_d, (unsigned long long)MAXL * BDIM));
	cudaErrchk(cudaMalloc(&T.new_edges_global_d, (unsigned long long)MAXL * EDGE_F * BDIM * sizeof(Edge)));
	cudaErrchk(cudaMalloc(&T.new_edge_bounds_global_d, (unsigned long long)(MAXL+1) * BDIM * sizeof(int)));
	cudaErrchk(cudaMalloc(&T.new_end_nodes_global_d, (unsigned long long)MAXL * BDIM));
	cudaErrchk(cudaMalloc(&T.new_sequence_ids_global_d, (unsigned long long)WL * MAXL * BDIM));

	cudaErrchk(cudaMalloc(&T.dyn_letters_global_d, (unsigned long long)MAXL * BDIM));
	cudaErrchk(cudaMalloc(&T.dyn_edges_global_d, (unsigned long long)MAXL * EDGE_F * BDIM * sizeof(Edge)));
	cudaErrchk(cudaMalloc(&T.dyn_edge_bounds_global_d, (unsigned long long)(MAXL+1) * BDIM * sizeof(int)));
	cudaErrchk(cudaMalloc(&T.dyn_end_nodes_global_d, (unsigned long long)MAXL * BDIM));
	cudaErrchk(cudaMalloc(&T.dyn_sequence_ids_global_d, (unsigned long long)WL * MAXL * BDIM));

	cudaErrchk(cudaMalloc(&T.moves_global_d, (unsigned long long)2 * (MAXL+1) * (SL+1) * BDIM * sizeof(unsigned char)));
	cudaErrchk(cudaMalloc(&T.diagonals_global_sc_d, (unsigned long long)(MAXL+1)*(SL+1) * BDIM * sizeof(short)));
	cudaErrchk(cudaMalloc(&T.diagonals_global_gx_d, (unsigned long long)(MAXL+1)*(SL+1) * BDIM * sizeof(short)));
	cudaErrchk(cudaMalloc(&T.diagonals_global_gy_d, (unsigned long long)(MAXL+1)*(SL+1) * BDIM * sizeof(short)));
	cudaErrchk(cudaMalloc(&T.d_offsets_global_d, (unsigned long long)(MAXL + SL+1) * BDIM * sizeof(int)));
	cudaErrchk(cudaMalloc(&T.x_to_ys_d, (unsigned long long)MAXL * BDIM * sizeof(int)));
	cudaErrchk(cudaMalloc(&T.y_to_xs_d, (unsigned long long)MAXL * BDIM * sizeof(int)));
	//cout << "Alloc completed\n";
}

inline void gpu_POA_free(TaskRefs &T, TaskType TTy){

	//cout << "Freeing memory for " << TTy << " \n";
/*
	cudaErrchk(cudaFree(T.sequences_d));
	cudaErrchk(cudaFree(T.seq_offsets_d));
	cudaErrchk(cudaFree(T.nseq_offsets_d));
	cudaErrchk(cudaFree(T.result_d));
	cudaErrchk(cudaFree(T.res_size_d));

	cudaErrchk(cudaFree(T.lpo_edge_offsets_d));
	cudaErrchk(cudaFree(T.lpo_letters_d));
	cudaErrchk(cudaFree(T.lpo_edges_d));
	cudaErrchk(cudaFree(T.edge_bounds_d));
	cudaErrchk(cudaFree(T.end_nodes_d));
	cudaErrchk(cudaFree(T.sequence_ids_d));

	cudaErrchk(cudaFree(T.new_letters_global_d));
	cudaErrchk(cudaFree(T.new_edges_global_d));
	cudaErrchk(cudaFree(T.new_edge_bounds_global_d));
	cudaErrchk(cudaFree(T.new_end_nodes_global_d));
	cudaErrchk(cudaFree(T.new_sequence_ids_global_d));
	
	cudaErrchk(cudaFree(T.dyn_letters_global_d));
	cudaErrchk(cudaFree(T.dyn_edges_global_d));
	cudaErrchk(cudaFree(T.dyn_edge_bounds_global_d));
	cudaErrchk(cudaFree(T.dyn_end_nodes_global_d));
	cudaErrchk(cudaFree(T.dyn_sequence_ids_global_d));

	cudaErrchk(cudaFree(T.moves_global_d));
	cudaErrchk(cudaFree(T.diagonals_global_sc_d));
	cudaErrchk(cudaFree(T.diagonals_global_gx_d));
	cudaErrchk(cudaFree(T.diagonals_global_gy_d));
	cudaErrchk(cudaFree(T.d_offsets_global_d));
	cudaErrchk(cudaFree(T.x_to_ys_d));
	cudaErrchk(cudaFree(T.y_to_xs_d));

	cudaErrchk(cudaFree(T.old_len_global_d));
	cudaErrchk(cudaFree(T.dyn_len_global_d));
*/
	cudaDeviceReset();
	
	free(T.result);
	free(T.res_size);

}

template<int SL, int MAXL, int WL, int BDIM> 
void gpu_POA(vector<Task<vector<string>>> &input, TaskRefs &T, vector<Task<vector<string>>> &result_GPU, int res_gpu_offs){
	
	int input_size = input.size();
	int N_BL = (input_size - 1) / BDIM + 1;
	int LAST_BATCH_SIZE = (input_size - 1) % BDIM + 1;
	
	vector<vector<string>> result_data(input_size);

	cout << "Executing POA: BDIM = " << BDIM << ", BATCHES = " << N_BL << endl;
	//cudaSetDevice(0);

	T.nseq_offsets = vector<int>(BDIM);

	cout << "Assign device mem\n";

	assign_device_memory<SL, MAXL, WL><<<1, 1>>>(T.lpo_edge_offsets_d, T.lpo_letters_d, T.lpo_edges_d, 
				       T.edge_bounds_d, T.end_nodes_d, T.sequence_ids_d, 
				       T.new_letters_global_d, T.new_edges_global_d, T.new_edge_bounds_global_d, 
				       T.new_end_nodes_global_d, T.new_sequence_ids_global_d, 
				       T.dyn_letters_global_d, T.dyn_edges_global_d, T.dyn_edge_bounds_global_d, 
				       T.dyn_end_nodes_global_d, T.dyn_sequence_ids_global_d,
				       T.moves_global_d, T.diagonals_global_sc_d, T.diagonals_global_gx_d, T.diagonals_global_gy_d, 
                                       T.d_offsets_global_d, T.x_to_ys_d, T.y_to_xs_d, T.old_len_global_d, T.dyn_len_global_d, BDIM);

	cudaStreamSynchronize(0);

	for(int b = 0; b < N_BL; b++){

		cout << "Batch: " << b << "\n";
		
		int block_offset = b * BDIM;
		int BLOCKS;
		if(b == N_BL-1){
			BLOCKS = LAST_BATCH_SIZE;
		}else{
			BLOCKS = BDIM;
		}

		//cout << "Init kernel parameters\n";
		init_kernel_block_parameters<BDIM>(input, &T.sequences, T.nseq_offsets, T.seq_offsets, &T.tot_nseq, block_offset);
		
		//cout << "Start memcpy\n";
		//cout << "Memcpy of " << T.seq_offsets[T.tot_nseq-1] << " bytes\n";
		cudaErrchk(cudaMemcpy(T.sequences_d, T.sequences, T.seq_offsets[T.tot_nseq-1], cudaMemcpyHostToDevice));
		cudaErrchk(cudaMemcpy(T.seq_offsets_d, T.seq_offsets.data(), T.tot_nseq * sizeof(int), cudaMemcpyHostToDevice));
		cudaErrchk(cudaMemcpy(T.nseq_offsets_d, T.nseq_offsets.data(), (unsigned long long)BLOCKS * sizeof(int), cudaMemcpyHostToDevice));
		
		//cout << "Compute edge offsets\n";

		compute_edge_offsets<SL, MAXL, WL><<<BLOCKS, WL>>>(T.seq_offsets_d, T.nseq_offsets_d);
		
		cudaStreamSynchronize(0);
		
		//cout << "Generate LPO\n";

		for(int i = 0; i < WL; i++){
			generate_lpo<SL, MAXL, WL><<<BLOCKS, SL+1>>>(T.sequences_d, T.seq_offsets_d, T.nseq_offsets_d, i);
		}
		
		int i_seq_idx = 0;
		for(int j_seq_idx = 1; j_seq_idx < WL; j_seq_idx++){
			
			//cout << "Alignment of sequence " << j_seq_idx << "\n";

			cudaStreamSynchronize(0);
				
			compute_d_offsets<SL, MAXL, WL><<<BLOCKS, 1>>>(i_seq_idx, j_seq_idx, T.nseq_offsets_d);
			
			cudaStreamSynchronize(0); 
			
			init_diagonals<SL, MAXL, WL><<<BLOCKS,1>>>(i_seq_idx, j_seq_idx, T.max_gapl, T.uses_global, T.nseq_offsets_d);
			
			//cout << "Alignment kernel call\n";

			sw_align<SL, MAXL, WL><<<BLOCKS, SL+1>>>(i_seq_idx, j_seq_idx, T.max_gapl, T.uses_global, T.nseq_offsets_d);
			
			cudaStreamSynchronize(0); 

			compute_new_lpo_size<SL, MAXL, WL><<<BLOCKS, 1>>>(i_seq_idx, j_seq_idx, T.nseq_offsets_d);
			
			cudaStreamSynchronize(0); 
			
			//cout << "Fusion kernel call\n";

			fuse_lpo<SL, MAXL, WL><<<BLOCKS, 1>>>(i_seq_idx, j_seq_idx, T.nseq_offsets_d);		
			
			cudaStreamSynchronize(0); 
			
			copy_new_lpo_data<SL, MAXL, WL><<<BLOCKS, MAXL+1>>>(j_seq_idx, T.nseq_offsets_d);
		
			//cout << "Copy of new LPO data inititated\n";
		}
		
		//cout << "Copy result sizes\n";

		copy_result_sizes<SL, MAXL, WL><<<BLOCKS, 1>>>(T.nseq_offsets_d, T.res_size_d);
		
		cudaStreamSynchronize(0); 
		
		suffix_sum<SL, MAXL, WL><<< 1, 1 >>>(T.res_size_d, BLOCKS);
		
		cudaStreamSynchronize(0); 
		
		//cout << "Computing final result\n";

		for(int i = 0; i < WL; i++){
			compute_result<SL, MAXL, WL><<<BLOCKS, MAXL>>>(T.nseq_offsets_d, T.result_d, T.res_size_d, i);
		}

		cudaStreamSynchronize(0);
		
		//cout << "Initiating memcpy back\n";

		cudaErrchk(cudaMemcpy(T.res_size, T.res_size_d, BDIM * sizeof(int), cudaMemcpyDeviceToHost));
		//cout << "mcpy2\n";
		cudaErrchk(cudaMemcpy(T.result, T.result_d, T.res_size[BLOCKS-1], cudaMemcpyDeviceToHost));
		
		//cout << "Window formation\n";

		result_GPU[input[block_offset].task_id].task_data = 
			form_window<MAXL>(T.result, T.nseq_offsets[0], T.res_size[0]/T.nseq_offsets[0]);
		result_GPU[input[block_offset].task_id].task_id = input[block_offset].task_id;
		result_GPU[input[block_offset].task_id].task_index = input[block_offset].task_index;

		int i = 1;
		auto i_it_end = input.begin() + block_offset + BLOCKS;
		for(auto it = input.begin() + block_offset+1; it != i_it_end; it++){

			int nseq_b = T.nseq_offsets[i] - T.nseq_offsets[i-1];
			int res_offset = T.res_size[i-1];
			int res_size_b = T.res_size[i] - T.res_size[i-1];
			char* res_offs = T.result + res_offset;
		
			result_GPU[it->task_id].task_id = it->task_id;
			result_GPU[it->task_id].task_index = it->task_index;
			result_GPU[it->task_id].task_data = form_window<MAXL>(res_offs, nseq_b, res_size_b/nseq_b);
			i++;
		}	

		free(T.sequences);	
	}
	
	cout << "Alignment task completed\n";	

	/*auto duration_alloc = duration_cast<microseconds>(alloc_end - alloc_begin);
	auto duration_exec = duration_cast<microseconds>(exec_end - alloc_end);
	auto duration_free = duration_cast<microseconds>(free_end - exec_end);
	auto duration_all = duration_cast<microseconds>(free_end - alloc_begin);

	cout << "alloc time = " << duration_alloc.count() << " microseconds" << endl;
	cout << "exec time = " << duration_exec.count() << " microseconds" << endl;
	cout << "free time = " << duration_free.count() << " microseconds" << endl;
	cout << "Total time = " << duration_all.count() << " microseconds" << endl;*/
}

}// end poa_gpu_utils

#endif //POA_GPU_H
