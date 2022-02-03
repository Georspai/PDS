#define EVALUATION_MODE 0

void start_state(int8_t** lat, int n);

int8_t** lattice_init(int n);

void read_bin(int8_t** lat, int n);

void resultCheck(int8_t**);

cudaError_t latticeInitCuda(int8_t*** lat, int8_t** dev_lat, int8_t** dev_lat_trans, int* n, int* k, int** dev_n, int** dev_k);

void printLattice(int8_t** lat, int* size);