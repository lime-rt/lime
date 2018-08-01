#ifndef LOCAL_ERR_H
#define LOCAL_ERR_H


#define ERR_STR_LEN	255

typedef struct  {
  int status;
  char message[ERR_STR_LEN+1];
} errType;

errType	init_local_err(void);
errType	write_local_err(int status, char *message);
void	error(int exitStatus, char *message);

#endif /* LOCAL_ERR_H */

