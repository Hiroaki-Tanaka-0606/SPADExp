// calcPSF log

#include <cstring>
#include <iostream>

#include "variables_ext.hpp"

using namespace std;
void write_log(char* message){
	strncpy(Log_buffer, message, Log_length+1);
	if(strlen(message)>Log_length){
		Log_buffer[Log_length]='\0';
	}
	if(Log_file_set){
	  fprintf(Log_file_obj, "%s\n",Log_buffer);
	}
	if(Console_log){
		cout << Log_buffer << endl;
	}
}
