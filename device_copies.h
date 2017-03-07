#ifndef DEVICE_COPIES_H
#define DEVICE_COPIES_H

/**
 * \brief function to do a host-to-device copy
 *
 * @param[in] dest 		- cuda pointer on device
 * @param[in] source 	- pointer on host
 * @param[in] bytes 	- number of bytes to copy
 */
void copy_to_device(void*,void*,unsigned);
/**
 * \brief function to do a device-to-host copy
 *
 * @param[in] dest 		- pointer on host
 * @param[in] source 	- cuda pointer on device
 * @param[in] bytes 	- number of bytes to copy
 */
void copy_from_device(void*,void*,unsigned);
/**
 * \brief function to do a device memory allocation
 *
 * @param[inout] dest 		- pointer on host
 * @param[inout] bytes 	- number of bytes to copy
 */
void allocate_on_device(void**,unsigned);
/**
 * \brief function to do a device memory allocation
 *
 * @param[inout] dest 		- pointer on device to deallocate
 */
void deallocate_on_device(void*);

#endif
