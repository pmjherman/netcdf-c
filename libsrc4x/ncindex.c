/*
Copyright (c) 1998-2017 University Corporation for Atmospheric Research/Unidata
See LICENSE.txt for license information.
*/

/** \file \internal
Internal netcdf-4 functions.

This file contains functions for manipulating ncindex objects.

Warning: This code depends critically on the assumption that
|void*| == |uintptr_t|

*/


/* Define this for debug so that table sizes are small */
#define SMALLTABLE

#include "config.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#ifdef HAVE_STDINT_H
#include <stdint.h>
#endif
#include <assert.h>

#include "includes4.h"

/* Keep the table sizes small initially */
#define DFALTTABLESIZE 7

#ifdef NCUSEHASH
extern void printhashmap(NC_hashmap*);
#endif

/* Locate object by name in an NCindex */
NC_OBJ*
ncindexlookup(NCindex* ncindex, const char* name)
{
   NC_OBJ* obj = NULL;
   if(ncindex == NULL || name == NULL)
	return NULL;
    {
#ifdef NCIUSEHASH
        uintptr_t index;
        assert(ncindex->map != NULL);
        if(!NC_hashmapget(ncindex->map,(void*)name,strlen(name),&index))
	    return NULL; /* not present */
	obj = (NC_OBJ*)nclistget(ncindex->list,(size_t)index);
#else
	int i;
        for(i=0;i<nclistlength(ncindex->list);i++) {
	    NC_OBJ* o = (NC_OBJ*)ncindex->list->content[i];
	    if(strcmp(o->name,name)==0) return o;
	}
#endif
    }
    return obj;
}

/* Get ith object in the vector */
NC_OBJ*
ncindexith(NCindex* index, size_t i)
{
   if(index == NULL) return NULL;
   assert(index->list != NULL);
   return nclistget(index->list,i);
}

/* See if x is contained in the index */
/* Return vector position if in index, otherwise return -1. */
int
ncindexfind(NCindex* index, NC_OBJ* nco)
{
    int i;
    NClist* list;
    if(index == NULL || nco == NULL) return -1;
    list = index->list;
    for(i=0;i<nclistlength(list);i++) {
	NC_OBJ* o = (NC_OBJ*)list->content[i];
	if(nco == o) return i;
    }
    return -1;
}

/* Add object to the end of the vector, also insert into the hashmaps */
/* Return 1 if ok, 0 otherwise.*/
int
ncindexadd(NCindex* ncindex, NC_OBJ* obj)
{
   if(ncindex == NULL) return 0;
   if(!nclistpush(ncindex->list,obj)) return 0;
#ifdef NCIUSEHASH
   {
   uintptr_t index; /*Note not the global id */
   index = (uintptr_t)nclistlength(ncindex->list);
   NC_hashmapadd(ncindex->map,index,(void*)obj->name,strlen(obj->name));
   }
#endif
   return 1;
}

/* Remove object from index by position in the vector.
 * WARNING: This will compress the vector by one, so
 * object id's may be affected. 
 * Return 1 if ok, 0 otherwise.*/
int
ncindexidel(NCindex* index, size_t i)
{
   NC_OBJ* obj;
   if(index == NULL) return 0;
   obj = nclistremove(index->list,i);
   if(obj == NULL)
	return 0; /* not present */
#ifdef NCIUSEHASH
   /* Remove from the hash map by deactivating its entry */
   if(!NC_hashmapdeactivate(index->map,(uintptr_t)i))
	return 0; /* not present */
#endif
   return 1;
}

/*Return a duplicate of the index's vector */
/* Return list if ok, NULL otherwise.*/
NC_OBJ**
ncindexdup(NCindex* index)
{
    if(index == NULL || nclistlength(index->list) == 0)
	return NULL;
    return (NC_OBJ**)nclistdup(index->list);
}

/*
Rebuild the list map by rehashing all entries
using their current, possibly changed id and name
*/
/* Return 1 if ok, 0 otherwise.*/
int
ncindexrehash(NCindex* index)
{ 
#ifdef NCIUSEHASH
    size_t i;
    size_t size = nclistlength(index->list);
    NC_OBJ** contents = (NC_OBJ**)nclistextract(index->list);
    /* Reset the list map */
    if(!ncindexclear(index))
	return 0;
    if(!ncindexinit(index,size))
	return 0;
    /* Now, reinsert all the attributes except NULLs */
    for(i=0;i<size;i++) {
	NC_OBJ* tmp = contents[i];
	if(tmp == NULL) continue; /* ignore */
	if(!ncindexadd(index,tmp))
	    return 0;
    }
#endif
    return 1;    
}

/* Free a list map */
int
ncindexfree(NCindex* index)
{
    if(index == NULL) return 1;
    nclistfree(index->list);
    NC_hashmapfree(index->map);
    free(index);    
    return 1;
}

/* Create a new index */
NCindex*
ncindexnew(size_t size0)
{
    NCindex* index = NULL;
    size_t size = (size0 == 0 ? DFALTTABLESIZE : size0);
    index = calloc(1,sizeof(NCindex));
    if(index == NULL) return NULL;
    index->list = nclistnew();
    if(index->list == NULL) {ncindexfree(index); return NULL;}
    nclistsetalloc(index->list,size);
#ifdef NCIUSEHASH
    index->map = NC_hashmapnew(size);
    if(index->map == NULL) {ncindexfree(index); return NULL;}
#endif
    return index;
}

#ifdef NCIUSEHASH
/* Debug: handle the data key part */
static const char*
keystr(NC_hentry* e)
{
    if(e->keysize < sizeof(uintptr_t))
	return (const char*)(&e->key);
    else 
	return (const char*)(e->key);
}
#endif

int
ncindexverify(NCindex* lm, int dump)
{
    size_t i;
    NClist* l = lm->list;
    int nerrs = 0;
#ifdef NCIUSEHASH
    size_t m;
    NC_hashmap* map = lm->map;
#endif

    if(lm == NULL) {
	fprintf(stderr,"index: <empty>\n");
	return 1;
    }  
    if(dump) {
	fprintf(stderr,"-------------------------\n");
#ifdef NCIUSEHASH
        if(lm->map->active == 0) {
	    fprintf(stderr,"hash: <empty>\n");
	    goto next1;
	}
	for(i=0;i < lm->map->alloc; i++) {
	    NC_hentry* e = &lm->map->table[i];
	    if(e->flags != 1) continue;
	    fprintf(stderr,"hash: %ld: data=%lu key=%s\n",(unsigned long)i,(unsigned long)e->data,keystr(e));
	    fflush(stderr);
	}
next1:
#endif
        if(nclistlength(l) == 0) {
	    fprintf(stderr,"list: <empty>\n");
	    goto next2;
	}
	for(i=0;i < nclistlength(l); i++) {
	    const char** a = (const char**)nclistget(l,i);
	    fprintf(stderr,"list: %ld: name=%s\n",(unsigned long)i,*a);
	    fflush(stderr);
	}
	fprintf(stderr,"-------------------------\n");
	fflush(stderr);
    }

next2:
#ifdef NCIUSEHASH
    /* Need to verify that every entry in map is also in vector and vice-versa */

    /* Verify that map entry points to same-named entry in vector */
    for(m=0;m < lm->map->alloc; m++) {
	NC_hentry* e = &lm->map->table[m];
        char** object = NULL;
	char* oname = NULL;
	uintptr_t udata = (uintptr_t)e->data;
	if((e->flags & 1) == 0) continue;
	object = nclistget(l,(size_t)udata);
        if(object == NULL) {
	    fprintf(stderr,"bad data: %d: %lu\n",(int)m,(unsigned long)udata);
	    nerrs++;
	} else {
	    oname = *object;
	    if(strcmp(oname,keystr(e)) != 0)  {
	        fprintf(stderr,"name mismatch: %d: %lu: hash=%s list=%s\n",
			(int)m,(unsigned long)udata,keystr(e),oname);
	        nerrs++;
	    }
	}
    }
    /* Walk vector and mark corresponding hash entry*/
    if(nclistlength(l) == 0 || lm->map->active == 0)
	goto done; /* cannot verify */
    for(i=0;i < nclistlength(l); i++) {
	int match;
	const char** xp = (const char**)nclistget(l,i);
        /* Walk map looking for *xp */
	for(match=0,m=0;m < lm->map->active; m++) {
	    NC_hentry* e = &lm->map->table[m];
	    if((e->flags & 1) == 0) continue;
	    if(strcmp(keystr(e),*xp)==0) {
		if((e->flags & 128) == 128) {
		    fprintf(stderr,"%ld: %s already in map at %ld\n",(unsigned long)i,keystr(e),(unsigned long)m);
		    nerrs++;
		}
		match = 1;
		e->flags += 128;
	    }
	}
	if(!match) {
	    fprintf(stderr,"mismatch: %d: %s in vector, not in map\n",(int)i,*xp);
	    nerrs++;
	}
    }
    /* Verify that every element in map in in vector */
    for(m=0;m < lm->map->active; m++) {
	NC_hentry* e = &lm->map->table[m];
	if((e->flags & 1) == 0) continue;
	if((e->flags & 128) == 128) continue;
	/* We have a hash entry not in the vector */
	fprintf(stderr,"mismatch: %d: %s->%lu in hash, not in vector\n",(int)m,keystr(e),(unsigned long)e->data);
	nerrs++;
    }
    /* clear the 'touched' flag */
    for(m=0;m < lm->map->active; m++) {
	NC_hentry* e = &lm->map->table[m];
	e->flags &= ~128;
    }

done:
#endif /*NCIUSEHASH*/
    fflush(stderr);
    return (nerrs > 0 ? 0: 1);
}

static const char*
sortname(NC_SORT sort)
{
    switch(sort) {
    case NCNAT: return "NCNAT";
    case NCVAR: return "NCVAR";
    case NCDIM: return "NCDIM";
    case NCATT: return "NCATT";
    case NCTYP: return "NCTYP";
    case NCGRP: return "NCGRP";
    default: break;
    }
    return "unknown";
}

void
printindexlist(NCindex* lm)
{
    int i;
    if(lm == NULL) {
	fprintf(stderr,"<empty>\n");
	return;
    }
    for(i=0;i<nclistlength(lm->list);i++) {
	NC_OBJ* o = (NC_OBJ*)nclistget(lm->list,i);
        fprintf(stderr,"[%ld] sort=%s name=|%s| id=%lu\n",
		(unsigned long)i,sortname(o->sort),o->name,(unsigned long)o->id);
    }
}

#ifdef NCIUSEHASH
void
printindexmap(NCindex* lm)
{
    if(lm == NULL) {
	fprintf(stderr,"<empty>\n");
	return;
    }
    printhashmap(lm->map);
}
#endif
