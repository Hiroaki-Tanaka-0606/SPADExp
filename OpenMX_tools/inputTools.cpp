// calcPSF tools for OpenMX
// inputTools
// see also OpenMX/source/Inputtools.c

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>

#include "inputTools.hpp"

using namespace std;

int bufSize=4096;
int keySize=128;
int valSize=128;
char** input_char;
int line_count;

char* remove_space(char* pointer){
	while(pointer[0]==' ' || pointer[0]=='\t'){
		if(pointer[0]=='\0'){
			return pointer;
		}
		pointer++;
	}
	return pointer;
}


void load_file(ifstream* fp){
	string line;
	int i;

	// count the number of lines
	line_count=0;
	while(getline(*fp, line)){
		// cout << line << endl;
		line_count++;
	}

	printf("The input file has %d rows.\n", line_count);

	input_char=new char*[line_count];
	for(i=0; i<line_count; i++){
		input_char[i]=new char[bufSize+1];
	}

	// go back to the top
	fp->clear();
	fp->seekg(0);
	int line_number=0;
	int actual_line_length;
	while(getline(*fp, line)){
		actual_line_length=line.length();
		if(actual_line_length>bufSize){
			printf("Warning: line %d is longer than %d characters\n", line_number+1, bufSize);
			actual_line_length=bufSize;
		}
		line.copy(input_char[line_number], actual_line_length);
		input_char[line_number][actual_line_length]='\0';
		// this process is necessary for loading vector values
		input_char[line_number]=remove_space(input_char[line_number]);
		line_number++;
	}
}

int load_logical(const char* key, bool* value){
	char keyBuf[keySize];
	char valBuf[valSize];
	int i, j;
	bool exists=false;
	const char* trueStrs[5];
	const char* falseStrs[5];
	int strChoices=5;
	trueStrs[0]="on";   trueStrs[1]="yes"; trueStrs[2]="true";   trueStrs[3]=".true.";   trueStrs[4]="ok";
	falseStrs[0]="off"; falseStrs[1]="no"; falseStrs[2]="false"; falseStrs[3]=".false."; falseStrs[4]="ng"; 
																					 
	
	for(i=0; i<line_count; i++){
		if(sscanf(input_char[i], "%s %s", keyBuf, valBuf)==2 && strcmp(keyBuf, key)==0){
			if(exists){
				printf("Error: key %s already exists\n", key);
				return 0;
			}
			exists=true;
			for(j=0; j<strChoices; j++){
				if(strncasecmp(valBuf, trueStrs[j], strlen(trueStrs[j]))==0){
					*value=true;
				}
				if(strncasecmp(valBuf, falseStrs[j], strlen(falseStrs[j]))==0){
					*value=false;
				}
			}
				
		}
	}
	if(exists){
		return 1;
	}else{
		printf("Error: key %s does not exist\n", key);
		return 0;
	}
	
}


int load_int(const char* key, int* value){
	char keyBuf[keySize];
	int valBuf;
	int i;
	bool exists=false;
	
	for(i=0; i<line_count; i++){
		if(sscanf(input_char[i], "%s %d", keyBuf, &valBuf)==2 && strcmp(keyBuf, key)==0){
			if(exists){
				printf("Error: key %s already exists\n", key);
				return 0;
			}
			exists=true;
			*value=valBuf;
		}
	}
	if(exists){
		return 1;
	}else{
		printf("Error: key %s does not exist or parse error happens\n", key);
		return 0;
	}
}

int load_doublev(const char* key, const int count, double* values){
	char keyBuf[keySize];
	char valBuf[valSize];
	int i, j;
	bool exists=false;
	
	for(i=0; i<line_count; i++){
		if(sscanf(input_char[i], "%s", keyBuf)==1 && strcmp(keyBuf, key)==0){
			if(exists){
				printf("Error: key %s already exists\n", key);
				return 0;
			}
			exists=true;
			char* current_pointer=input_char[i]+strlen(keyBuf);
			current_pointer=remove_space(current_pointer);
			for(j=0; j<count; j++){
				if(sscanf(current_pointer, "%lf", &values[j])!=1){
					printf("Error in parsing key %s\n", key);
					return 0;
				}else{
					sscanf(current_pointer, "%s", valBuf);
					current_pointer+=strlen(valBuf);
					current_pointer=remove_space(current_pointer);
				}
			}
		}
	}
	if(exists){
		return 1;
	}else{
		printf("Error: key %s does not exist\n", key);
		return 0;
	}
}

int find_str(const char* key){
	char keyBuf[keySize];
	int i;
	bool exists=false;
	int line_number;
	
	for(i=0; i<line_count; i++){
		if(sscanf(input_char[i], "%s", keyBuf)==1 && strcmp(keyBuf, key)==0){
			if(exists){
				printf("Error: key %s already exists\n", key);
				return -1;
			}
			exists=true;
			line_number=i;
		}
	}
	if(exists){
		return line_number;
	}else{
		printf("key %s does not exist\n", key);
		return -1;
	}
}

char* get_line(int line_number){
	return input_char[line_number];
}

void copyInput(ofstream* fp){
	// copy the input to fp
	// neglect MO.fileout, MO.Nkpoint, and <MO.kpoint> block

	char keyBuf[keySize];
	int i;
	int line_MOk_start=find_str("<MO.kpoint");
	int line_MOk_end=find_str("MO.kpoint>");
	for(i=0; i<line_count; i++){
		if(!(line_MOk_start <= i && i <= line_MOk_end)){
			if(sscanf(input_char[i], "%s", keyBuf)==1 &&
				 (strcmp(keyBuf, "MO.fileout")==0 || strcmp(keyBuf, "MO.Nkpoint")==0)){
				*fp << "# ";
			}
		}else{
			*fp << "# ";
		}
		*fp << input_char[i] << endl;
	}
}
