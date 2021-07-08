// calcPSF parse inputs

int parse_char(char* line, char* variable_name, char* variable, int variable_size, char* value_buffer);
int parse_bool(char* line, bool* variable, char* value_buffer);
int parse_double(char* line, double* variable);
int parse_int(char* line, int* variable);
