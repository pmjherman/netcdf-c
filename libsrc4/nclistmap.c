/*
Copyright (c) 1998-2017 University Corporation for Atmospheric Research/Unidata
See LICENSE.txt for license information.
*/

/** \file \internal
Internal netcdf-4 functions.

This file contains functions for manipulating NC_listmap
objects.

Warning: This code depends critically on the assumption that
|void*| == |uintptr_t|

*/

#include "config.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#ifdef HAVE_STDINT_H
#include <stdint.h>
#endif
#include <assert.h>
#include "nclistmap.h"
#include "nc4internal.h"

/* Keep the table sizes small initially */
#define DFALTTABLESIZE 7

extern void printhashmap(NC_hashmap*);

/* Locate object by name in an NC_LISTMAP */
NC_OBJ*
NC_listmap_get(NC_listmap* listmap, const char* name)
{
   uintptr_t index;
   NC_OBJ* obj;
   if(listmap == NULL || name == NULL)
	return NULL;
   assert(listmap->map != NULL);
   if(!NC_hashmapget(listmap->map,(void*)name,strlen(name),&index))
	return NULL; /* not present */
   obj = (NC_OBJ*)nclistget(listmap->list,(size_t)index);
   return obj;
}

#if 0
/* Locate object by id in an NC_LISTMAP; also can be used to test if object is in this listmap */
NC_OBJ*
NC_listmap_iget(NC_listmap* listmap, size_t id)
{
   uintptr_t index;
   NC_OBJ* obj;
   if(listmap == NULL)
	return NULL;
   assert(listmap->idmap != NULL);
   if(!NC_hashmapget(listmap->idmap,&id,sizeof(id),&index))
	return NULL; /* not present */
   obj = (NC_OBJ*)nclistget(listmap->list,(size_t)index);
   return obj;
}
#endif

/* Get ith object in the vector */
void*
NC_listmap_ith(NC_listmap* listmap, size_t index)
{
   if(listmap == NULL) return NULL;
   assert(listmap->list != NULL);
   return nclistget(listmap->list,index);
}

/* Add object to the end of the vector, also insert into the hashmaps */
/* Return 1 if ok, 0 otherwise.*/
int
NC_listmap_add(NC_listmap* listmap, NC_OBJ* obj)
{
   uintptr_t index; /*Note not the global id */
   if(listmap == NULL) return 0;
   index = (uintptr_t)nclistlength(listmap->list);
   if(!nclistpush(listmap->list,obj)) return 0;
   NC_hashmapadd(listmap->map,index,(void*)obj->name,strlen(obj->name));
   return 1;
}

/* Remove object from listmap by index in the vector.
 * WARNING: This will compress the vector by one, so
 * object id's may be affected. 
 * Return 1 if ok, 0 otherwise.*/
int
NC_listmap_idel(NC_listmap* listmap, size_t index)
{
   NC_OBJ* obj;
   if(listmap == NULL) return 0;
   obj = nclistremove(listmap->list,index);
   if(obj == NULL)
	return 0; /* not present */
   /* Remove from the hash map by deactivating its entry */
   if(!NC_hashmapdeactivate(listmap->map,(uintptr_t)index))
	return 0; /* not present */
   return 1;
}

#if 0
/* Add object at specific index; will overwrite anything already there;
   obj == NULL is ok.
   Return 1 if ok, 0 otherwise.
 */
int
NC_listmap_iput(NC_listmap* listmap, size_t pos, void* obj)
{
    if(listmap == NULL) return 0;
    if(obj != NULL) {
	uintptr_t data = 0;
	const char* name = *(const char**)obj;
	if(name != NULL) return 0;
	if(pos >= nclistlength(listmap->list)) return 0;
	/* Temporarily remove from hashmap */
	if(!NC_hashmapremove(listmap->map,(void*)name,strlen(name),NULL))
	    return 0; /* not there */
	/* Reinsert with new pos */
	data = (uintptr_t)pos;
	NC_hashmapadd(listmap->map,(void*)data,*(char**)obj);
    }
    /* Insert at pos into listmap vector */
    nclistset(listmap->list,pos,obj);
    return 1;
}
#endif

#if 0
/* Remove object from listmap; assume cast (char**)target is defined */
/* WARNING: This will leave a NULL hole in the vector */
/* Return 1 if ok, 0 otherwise.*/
int
NC_listmap_del(NC_listmap* listmap, NC_OBJ* target)
{
   uintptr_t pos;
   if(listmap == NULL || target == NULL) return 0;
   /* Remove from the hash maps */
   if(!NC_hashmapremove(listmap->idmap,&target->id,sizeof(target->id),&pos))
	return 0; /* not present */
   if(!NC_hashmapremove(listmap->map,target->name,strlen(target->name),NULL))
	return 0; /* not present */
   /* NULLify entry in the vector; we do this so we do not have to rehash
      all higher entries */
   nclistset(listmap->list,(size_t)pos,NULL);
   return 1;
}

/* Pseudo iterator; start index at 0, return 0 when complete.
   Usage:
      size_t iter;
      NC_OBJ** data	  
      for(iter=0;NC_listmap_next(listmap,iter,&data);iter++) {f(data);}
*/
size_t
NC_listmap_next(NC_listmap* listmap, size_t index, NC_OBJ** datap)
{
    size_t len = nclistlength(listmap->list);
    if(datap) *datap = NULL;
    if(len == 0) return 0;
    if(index >= nclistlength(listmap->list)) return 0;
    if(datap) *datap = (NC_OBJ*)nclistget(listmap->list,index);
    return index+1;
}

/* Reverse pseudo iterator; start index at 0, return 1 if more data, 0 if done.
   Differs from NC_listmap_next in that it iterates from last to first.
   This means that the iter value cannot be directly used as an index
   for e.g. NC_listmap_iget().
   Usage:
      size_t iter;
      void* data;
      for(iter=0;NC_listmap_next(listmap,iter,&data);iter++) {f(data);}
*/
size_t
NC_listmap_prev(NC_listmap* listmap, size_t iter, void** datap)
{
    size_t len = nclistlength(listmap->list);
    size_t index;
    if(datap) *datap = NULL;
    if(len == 0) return 0;
    if(iter >= len) return 0;
    index = (len - iter) - 1;
    if(datap) *datap = nclistget(listmap->list,index);
    return iter+1;
}
#endif

/*Return a duplicate of the listmap's vector */
/* Return list if ok, NULL otherwise.*/
NC_OBJ**
NC_listmap_dup(NC_listmap* listmap)
{
    if(listmap == NULL || nclistlength(listmap->list) == 0)
	return NULL;
    return (NC_OBJ**)nclistdup(listmap->list);
}

/*
Rebuild the list map by rehashing all entries
using their current, possibly changed id and name
*/
/* Return 1 if ok, 0 otherwise.*/
int
NC_listmap_rehash(NC_listmap* listmap)
{ 
    size_t i;
    size_t size = nclistlength(listmap->list);
    NC_OBJ** contents = (NC_OBJ**)nclistextract(listmap->list);
    /* Reset the list map */
    if(!NC_listmap_clear(listmap))
	return 0;
    if(!NC_listmap_init(listmap,size))
	return 0;
    /* Now, reinsert all the attributes except NULLs */
    for(i=0;i<size;i++) {
	NC_OBJ* tmp = contents[i];
	if(tmp == NULL) continue; /* ignore */
	if(!NC_listmap_add(listmap,tmp))
	    return 0;
    }
    return 1;    
}

/* Clear a list map without free'ing the map itself */
int
NC_listmap_clear(NC_listmap* listmap)
{
    nclistfree(listmap->list);
    NC_hashmapfree(listmap->map);
    listmap->list = NULL;    
    listmap->map = NULL;    
    return 1;
}

/* Initialize a list map without malloc'ing the map itself */
int
NC_listmap_init(NC_listmap* listmap, size_t size0)
{
    size_t size = (size0 == 0 ? DFALTTABLESIZE : size0);
    listmap->list = nclistnew();
    if(listmap->list != NULL)
	nclistsetalloc(listmap->list,size);
    listmap->map = NC_hashmapnew(size);
    return 1;
}

/* Handle the data key part */
static const char*
keystr(NC_hentry* e)
{
    if(e->keysize < sizeof(uintptr_t))
	return (const char*)(&e->key);
    else 
	return (const char*)(e->key);
}

int
NC_listmap_verify(NC_listmap* lm, int dump)
{
    size_t i;
    NC_hashmap* map = lm->map;
    NClist* l = lm->list;
    size_t m;
    int nerrs = 0;

    if(dump) {
	fprintf(stderr,"-------------------------\n");
        if(map->active == 0) {
	    fprintf(stderr,"hash: <empty>\n");
	    goto next1;
	}
	for(i=0;i < map->alloc; i++) {
	    NC_hentry* e = &map->table[i];
	    if(e->flags != 1) continue;
	    fprintf(stderr,"hash: %ld: data=%lu key=%s\n",(unsigned long)i,(unsigned long)e->data,keystr(e));
	    fflush(stderr);
	}
next1:
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
    /* Need to verify that every entry in map is also in vector and vice-versa */

    /* Verify that map entry points to same-named entry in vector */
    for(m=0;m < map->alloc; m++) {
	NC_hentry* e = &map->table[m];
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
    if(nclistlength(l) == 0 || map->active == 0)
	goto done; /* cannot verify */
    for(i=0;i < nclistlength(l); i++) {
	int match;
	const char** xp = (const char**)nclistget(l,i);
        /* Walk map looking for *xp */
	for(match=0,m=0;m < map->active; m++) {
	    NC_hentry* e = &map->table[m];
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
    for(m=0;m < map->active; m++) {
	NC_hentry* e = &map->table[m];
	if((e->flags & 1) == 0) continue;
	if((e->flags & 128) == 128) continue;
	/* We have a hash entry not in the vector */
	fprintf(stderr,"mismatch: %d: %s->%lu in hash, not in vector\n",(int)m,keystr(e),(unsigned long)e->data);
	nerrs++;
    }
    /* clear the 'touched' flag */
    for(m=0;m < map->active; m++) {
	NC_hentry* e = &map->table[m];
	e->flags &= ~128;
    }

done:
    fflush(stderr);
    return (nerrs > 0 ? 0: 1);
}

static const char*
sortname(NC_SORT sort)
{
    switch(sort) {
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
printlistmaplist(NC_listmap* lm)
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

void
printlistmapmap(NC_listmap* lm)
{
    if(lm == NULL) {
	fprintf(stderr,"<empty>\n");
	return;
    }
    printhashmap(lm->map);
}
