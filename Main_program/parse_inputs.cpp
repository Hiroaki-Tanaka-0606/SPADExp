// calcPSF parse inputs

#include <cstring>
#include <iostream>

// parse value of char*
int parse_char(char* line, char* variable_name, char* variable, int variable_size, char* value_buffer){
	int sscanf_status=sscanf(line, "%*s %s", value_buffer);
	if(sscanf_status<1){
		// parse failed
		return 0;
	}
	strncpy(variable, value_buffer, variable_size+1);
	if(strlen(value_buffer)>variable_size){
		printf("Warning: length of the value of %s is longer than %d, the last part is ignored\n", variable_name, variable_size);
		variable[variable_size]='\0';
	}
	return 1;
}

int parse_bool(char* line, bool* variable, char* value_buffer){
	int sscanf_status=sscanf(line, "%*s %s", value_buffer);
	if(sscanf_status<1){
		// parse failed
		return 0;
	}
	if(strcmp(value_buffer, "True")==0 || strcmp(value_buffer, "TRUE")==0 || strcmp(value_buffer, "true")==0){
		*variable=true;
		return 1;
	}
	if(strcmp(value_buffer, "False")==0 || strcmp(value_buffer, "FALSE")==0 || strcmp(value_buffer, "false")==0){
		*variable=false;
		return 1;
	}
	return 0;
}
int parse_double(char* line, double* variable){
	int sscanf_status=sscanf(line, "%*s %lf", variable);
	if(sscanf_status<1){
		// parse failed
		return 0;
	}
	return 1;
}
int parse_int(char* line, int* variable){
	int sscanf_status=sscanf(line, "%*s %d", variable);
	if(sscanf_status<1){
		// parse failed
		return 0;
	}
	return 1;
}

