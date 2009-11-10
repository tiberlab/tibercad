#ifndef _DATABASECALL_H_
#define _DATABASECALL_H_


#include <string>
#include <iostream>  
#include "sqlite3.h"
#include <stdio.h>
#include <sstream>

#include <vector>
#include <string>


#define MAX_TABLE_NAME_LENGTH 64
#define MAX_FIELD_NAME_LENGTH 64
#define MAX_PROPERTY_NAME_LENGTH 256
#define MAX_FIELD_VALUE_LENGTH 256
#define MAX_SQL_TEXT_LENGTH 1024
#define MAX_ID_LIST_LENGTH 512



using namespace   std;





class DataBaseCall
{
 public:


  DataBaseCall(const char* filename);
  //DataBaseCall::DataBaseCall(string material, string structure ,  string sim_class,  string model);



  ~DataBaseCall();

  void set_query( vector<string>& vector_query);

  // double db_class::get_data(string& prop);

  string get_data(const string  prop);



  private:

  const char*  db_location;

  vector<string> prop_label;
  vector<string> val1_value;


  char *zErrMsg; //= 0;
  char selectedStructures[MAX_ID_LIST_LENGTH];
  char selectedStructureName[MAX_ID_LIST_LENGTH];
  char selectedModels[MAX_ID_LIST_LENGTH];
  char selectedModels1[MAX_ID_LIST_LENGTH];
  char selectedModels2[MAX_ID_LIST_LENGTH];
  char selectedModels4[MAX_ID_LIST_LENGTH];
  char selectedClasses[MAX_ID_LIST_LENGTH];
  int rc;







  // Selected material
  char Material[MAX_FIELD_VALUE_LENGTH];
  // Selected material id read by callbackGetMaterialId
  char MaterialId[MAX_FIELD_VALUE_LENGTH];
  // Comma-separated list of ids of tables or properties selected by user
  char SelectedIds[MAX_ID_LIST_LENGTH];

  // Current database
  sqlite3 *Db;



  bool   exec(const char *sqlText, sqlite3_callback callbackFun, void *parm);

  // callback function for reading MaterialId
  static int callbackGetMaterialId(void *NotUsed, int argc, char **argv, char **azColName);

  // Callback function for further property printing
  static int callbackSelectPropertyValue(void *NotUsed, int argc, char **argv, char **azColName);


  static int callbackPrintDisplayName(void *print, int argc, char **argv, char **azColName);

  static int callbackGetStructureId(void *NotUsed, int argc, char **argv, char **azColName);

  bool GetStructure();

  bool GetMaterial();


  bool ExecuteSelect(
                     const char *sqlTemplate,
                     sqlite3_callback callbackFunction,
                     ...
                     );


  bool SelectIds0(
                  char *targetBuffer,
                  bool print,
                  const char *sqlTemplate, //format string for variable-length arguments list.
                  va_list marker
                  );


  bool SelectIds(
                 char *targetBuffer,
                 const char *sqlTemplate, //format string for variable-length arguments list.
                 ...
                 );




  bool UserSelectIds(
                     char *targetBuffer, 
                     const char *prompText, // prompt displayed to the user
                     const char *sqlTemplate, //format string for variable-length arguments list.
                     ...
                     );





  void query_db();



};


#endif /*_DATABASECALL_H_*/
