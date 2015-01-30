#ifndef DEVICE_COPIES_H
#define DEVICE_COPIES_H

void copy_to_device(void*,void*,unsigned);
void copy_from_device(void*,void*,unsigned);
void allocate_on_device(void**,unsigned);
void deallocate_on_device(void*);

#endif
