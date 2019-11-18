
struct record {
   
   double key;
   double col1;
   double col2;
   double col3;
   double col4;
   double col5;
   double col6;
   double col7;
   struct record* prev;
   struct record* next;
};

void append_node(struct record *lnode);
void insert_node(struct record *lnode, struct record *after);
void remove_node(struct record *lnode);
struct record* insert(struct record *p, struct record *list);

struct record *head, *tail;

int sortFile(char* filename) {
   
   FILE *fp;
   struct record *row, *curr;

   head = NULL;

   fp = fopen(filename, "r");   
   
   if(!fp) { 
      
      printf("\n\nUnable to open file to read - %s\n\n", filename);
      return -1;
   }
   
   while(!feof(fp)) {
      
      row = (struct record*) malloc(sizeof(struct record));
      
      row->prev = NULL;
      row->next = NULL;
      
      fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf\n", &row->col1, &row->col2, &row->col3, &row->col4, &row->col5, &row->col6, &row->col7);
      row->key = row->col1;
      
      head = insert(row, head);      
   }
   
   fclose(fp);
   
   fp = fopen(filename, "w");

   if(!fp) {
      
      printf("\n\nUnable to open file to write sorted data - %s\n\n", filename);
      return -1;
   }
   
   /* write the data back into the file */
   for(row = head; row != NULL; row = row->next) {
      
     fprintf(fp, "%e %e %e %e %e %e %e\n", row->col1, row->col2, row->col3, row->col4, row->col5, row->col6, row->col7);
   }
   
   fclose(fp);
   
   return 0;
}


struct record *insert(struct record *p, struct record *list) {
   
   struct record *q;
   
   /* first, we handle the case where data' should be the first element */
   if(list == NULL) {

      /* this is the first element */
      p->next = NULL;
      p->prev = NULL;
      
      return p;
      
   } else if (list->key > p->key) {
      
      /* now data should [be|becomes] the first element */
      p->next = list;
      p->prev = NULL;
      p->next->prev = p;
      
      return p;

   } else {
      
      /* search the linked list for the right location */
      q = list;
      while(q->key < p->key) {
	 
	 if(q->next != NULL) {
	    
	    q = q->next;
	    
	 } else {
	    
	    break;
	 }
      }
      
      if(q->next == NULL) {

	 if(q->key > p->key) {
	    
	    p->next = q;
	    p->prev = q->prev;
	    p->next->prev = p;
	    p->prev->next = p;

	 } else { 
	    
	    p->next = NULL;
	    p->prev = q;
	    p->prev->next = p;
	 }

      } else {
	 
	 p->prev = q->prev;
	 q->prev->next = p;
	 p->next = q;
	 p->next->prev = p;
      }
      
      return list;
   }
   
}

