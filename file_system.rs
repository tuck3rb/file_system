// NOTES FROM CLASS
// ~~~~~~~~~~~~~~~~

// Binary and
// &

// Binary or
// |

// Left shift
// 1 << 3
// Number << spaces moved left

// Right shift
// u16 >> 3
// Number >> spaces moved right


// Root directory:
// 	- Inode 0
// 	- Set bit 0 in inode block to 1

// read(block)
// write(block)

// u8::MAX     -   11111111
// U16& - bitwise and

// Create file
// 	- pick inode
// 	- Find lowest 0 bit
// 		- Inodeexist -> buffer
// 			- Examine each byte
// 				- < u8 .. MAX
// 				- Byte % 2
// 				- Byte / 2 	>> 1
				
				
				
				
// Offset into file
// 	- # bytes in file !
// 	- Add bytesoleafter array
// 		- Track total size!
// 	- Fill block
// 		- Write block to disk
// 		- Request now data block
// 		- Update inode
// 		- Clear bytes?

// MOST IMPORTANT
// ~~~~~~~~~~~~~~
// Blocksize (*8 is ceiling)
// Max#ofblocks
// Length of filename
// Max blocks per file --- --- --- --- --- --- [block 1, block 2, ...]  (2 + Max_Blocks_file = number of bytes in an inode)
// Max files stored
    // -> Relates to the total inodes to store
// Max file size < 2^16 bytes

// Copy test code into this (aka lib.rs? - do we rename this to lib.rs?)

// Block 0 : Inodes in use
// 0 - 000000001
// 1 - ?

// Block 1 : Datablocks in use
// 0 - 00111111
// 1 - 00000000

// File content buffer (load entire inode table into that buffer), then just do array read and writes directly into that table
// 

// Directory stuff
// Length: 8
// 'a', 'l', 'p', 'h', 'a', '\0', '\0', '\0'        (\o is null)
// 'b', 'e', 't', 'a', '\0', '\0', '\0', '\0'

// use bitwise or to set a given bit to be a 1
            // buffer[0]
            // want buffer to have a 1 in bit 0
            // buffer[0] 
            // buffer[0] |= 1
            // |   (bitwise or)

            // to set one
            // buffer[0] |= 2
            // set 2
            // buffer[0] |= 4

            // buffer[0] |= 1 << bit

            // to modify inode 19:
            // 19 /8 = 2 so we will be in 2 in the arrray of inodes or whatever
            // 19 % 2 = 3 so we will do the left shift to 3

            // buffer[2] |= 1 << 3


#![cfg_attr(not(test), no_std)]

#[derive(Copy, Clone, PartialEq, Eq)]
pub enum FileSystemResult<T: Copy + Clone> {
    Ok(T),
    Err(FileSystemError),
}

impl<T: Copy + Clone> FileSystemResult<T> {
    pub fn unwrap(&self) -> T {
        match self {
            FileSystemResult::Ok(v) => *v,
            FileSystemResult::Err(e) => panic!("Error: {e:?}"),
        }
    }
}

#[derive(Copy, Clone, Eq, PartialEq, Debug)]
pub enum FileSystemError {
    FileNotFound,
    FileNotOpen,
    NotOpenForRead,
    NotOpenForWrite,
    TooManyOpen,
    TooManyFiles,
    AlreadyOpen,
    DiskFull,
    FileTooBig,
    FilenameTooLong,
    InvalidFileDescriptor(usize),
}

#[derive(Debug, Copy, Clone)]
pub struct FileInfo<const MAX_BLOCKS: usize, const BLOCK_SIZE: usize> {
    inode: Inode<MAX_BLOCKS, BLOCK_SIZE>,
    inode_num: usize,
    current_block: usize,
    offset: usize,
    writing: bool,
    block_buffer: [u8; BLOCK_SIZE],
}

impl<const MAX_BLOCKS: usize, const BLOCK_SIZE: usize> FileInfo<MAX_BLOCKS,BLOCK_SIZE> {
    pub fn new_open_write(inode:Inode<MAX_BLOCKS, BLOCK_SIZE>, inode_num: usize) -> Self {
        Self {
            inode,
            inode_num,
            current_block: 0,
            offset: 0,
            writing: true,
            block_buffer: [0; BLOCK_SIZE],
        }
    }

    // pub fn new_open_read() -> Self {

    // }
}

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub struct Inode<const MAX_BLOCKS: usize, const BLOCK_SIZE: usize> {
    bytes_stored: u16,
    blocks: [u8; MAX_BLOCKS],
}

const INODE_FULL_BLOCK: usize = 0;
const DATA_FULL_BLOCK: usize = INODE_FULL_BLOCK + 1;
const INODE_TABLE_START: usize = DATA_FULL_BLOCK + 1;

#[derive(core::fmt::Debug)]
pub struct FileSystem<
    const MAX_OPEN: usize,
    const BLOCK_SIZE: usize,
    const NUM_BLOCKS: usize,
    const MAX_FILE_BLOCKS: usize,
    const MAX_FILE_BYTES: usize,
    const MAX_FILES_STORED: usize,
    const MAX_FILENAME_BYTES: usize,
> {
    open: [Option<FileInfo<MAX_FILE_BLOCKS, BLOCK_SIZE>>; MAX_OPEN],
    disk: ramdisk::RamDisk<BLOCK_SIZE, NUM_BLOCKS>,
    block_buffer: [u8; BLOCK_SIZE],
    file_content_buffer: [u8; MAX_FILE_BYTES],
    open_inodes: [bool; MAX_FILES_STORED],
}

impl<
        const MAX_OPEN: usize,
        const BLOCK_SIZE: usize,
        const NUM_BLOCKS: usize,
        const MAX_FILE_BLOCKS: usize,
        const MAX_FILE_BYTES: usize,
        const MAX_FILES_STORED: usize,
        const MAX_FILENAME_BYTES: usize,
    >
    FileSystem<
        MAX_OPEN,
        BLOCK_SIZE,
        NUM_BLOCKS,
        MAX_FILE_BLOCKS,
        MAX_FILE_BYTES,
        MAX_FILES_STORED,
        MAX_FILENAME_BYTES,
    >
{
    pub fn new(disk: ramdisk::RamDisk<BLOCK_SIZE, NUM_BLOCKS>) -> Self {
        assert_eq!(MAX_FILE_BYTES, MAX_FILE_BLOCKS * BLOCK_SIZE);
        assert!(NUM_BLOCKS <= u8::MAX as usize);
        assert!(MAX_FILE_BYTES <= u16::MAX as usize);
        let block_bits = BLOCK_SIZE * 8;
        assert!(MAX_FILES_STORED <= block_bits);
        assert!(MAX_FILES_STORED <= u16::MAX as usize);
        let result = Self {
            open: [None; MAX_OPEN],
            disk,
            block_buffer: [0; BLOCK_SIZE],
            file_content_buffer: [0; MAX_FILE_BYTES],
            open_inodes: [false; MAX_FILES_STORED],
        };
        assert!(result.num_inode_blocks() * 2 < NUM_BLOCKS);
        assert!(result.num_data_blocks() <= block_bits);
        assert_eq!(
            result.num_data_blocks() + result.num_inode_blocks() + 2,
            NUM_BLOCKS
        );
        assert!(result.num_inode_entries() <= u16::MAX as usize);
        assert!(result.num_inode_blocks() <= MAX_FILE_BLOCKS);
        result
    }

    pub fn assert_block(&self, block: usize, offset: usize, block_segment: &[u8]) {
        assert!(block < self.disk.num_blocks());
        let mut bytes = [0; BLOCK_SIZE];
        self.disk.read(block, &mut bytes);
        assert_eq!(block_segment, &bytes[offset..offset + block_segment.len()]);
    }

    pub fn max_file_size(&self) -> usize {
        MAX_FILE_BLOCKS * BLOCK_SIZE
    }

    pub fn num_inode_bytes(&self) -> usize {
        2 + MAX_FILE_BLOCKS
    }

    pub fn inodes_per_block(&self) -> usize {
        BLOCK_SIZE / self.num_inode_bytes()
    }

    pub fn num_inode_blocks(&self) -> usize {
        MAX_FILES_STORED / self.inodes_per_block()
    }

    pub fn num_data_blocks(&self) -> usize {
        NUM_BLOCKS - self.num_inode_blocks() - 2
    }

    pub fn num_inode_entries(&self) -> usize {
        self.inodes_per_block() * self.num_inode_blocks() * self.num_inode_bytes()
    }

    pub fn first_data_block(&self) -> usize {
        2 + self.num_inode_blocks()
    }

    pub fn open_read(&mut self, filename: &str) -> FileSystemResult<usize> {
        // Load the directory file and find the file's inode
        let file_inode_num = match self.find_file_inode(filename) {
            Some(inode_num) => inode_num,
            None => return FileSystemResult::Err(FileSystemError::FileNotFound),
        };
    
        // Check if inode is already open
        if self.open_inodes[file_inode_num] {
            return FileSystemResult::Err(FileSystemError::AlreadyOpen);
        }
    
        // Get inode from inode table
        let inode = self.get_inode_from_table(file_inode_num);
    
        // Create a file table entry for the newly opened file
        let mut file_info = FileInfo {
            inode: inode,
            inode_num: file_inode_num,
            current_block: 0,
            offset: 0,
            writing: true,
            block_buffer: [0; BLOCK_SIZE],
        };
    
        // Read in the first block of the newly opened file into the file table entry's buffer
        if file_info.inode.blocks[0] != 0 {
            self.disk.read(file_info.inode.blocks[0] as usize, &mut file_info.block_buffer);
        }
    
        // Mark the inode as open in open_inodes
        self.open_inodes[file_inode_num] = true;
    
        // Add the FileInfo to the open table and return the file descriptor
        let file_descriptor = self.add_to_open_table(file_info).unwrap(); 

        // Return the file descriptor
        FileSystemResult::Ok(file_descriptor)
        // todo!();

    }
    

    pub fn open_create(&mut self, filename: &str) -> FileSystemResult<usize> {
        let mut inode_block_buffer = [0; BLOCK_SIZE];
        let mut inode_table = [0; MAX_FILE_BYTES];
        let mut data_blocks_in_use = [0; BLOCK_SIZE];
        let mut directory_contents = [0; MAX_FILE_BYTES]; 
    
        // Read the inode and blocks from the disk
        self.disk.read(0, &mut inode_block_buffer);
        self.disk.read(1, &mut data_blocks_in_use);

        let iterations = (inode_table.len() + self.block_buffer.len() - 1) / self.block_buffer.len();
        for i in 0..iterations {  
            self.disk.read(0, &mut self.block_buffer); 
            for j in 0..self.block_buffer.len() {
                let index = i * self.block_buffer.len() + j;
                if index < inode_table.len() {
                    inode_table[index] = self.block_buffer[j];
                } else {
                    break;
                }
            }
        }

        let mut file_descriptor: usize = 0;
        let mut file_table_entry: FileInfo<MAX_FILE_BLOCKS, BLOCK_SIZE>;

        // Create the directory file if it does not exist
        if inode_block_buffer[0] == 0 {
            // Set the bit for inode 0 to 1
            inode_block_buffer[0] = 1;
    
            // Select the first data block for the directory
            for i in 0..=(self.first_data_block()) {
                data_blocks_in_use[i / 8] |= 1 << (i % 8); 
            }
    
            // Create an inode for the directory
            let mut directory_inode = Inode {
                bytes_stored: (MAX_FILENAME_BYTES * 2) as u16,
                blocks: [0; MAX_FILE_BLOCKS],
            };
            directory_inode.blocks[0] = self.first_data_block() as u8;
            // Save it in the inode table
            self.write_inode_to_table(0, &directory_inode);
            
            let new_inode_num = Self::get_inode(&mut inode_block_buffer);
            let data_block_num = Self::get_inode(&mut data_blocks_in_use);

            // create inode for file
            let mut new_inode = Inode::<MAX_FILE_BLOCKS, BLOCK_SIZE> {
                bytes_stored: 0,
                blocks: [0; MAX_FILE_BLOCKS],
            };
            new_inode.blocks[0] = data_block_num as u8;


            let file_table_entry = FileInfo {
                inode: new_inode,
                inode_num: new_inode_num,
                current_block: 0,
                offset: 0,
                writing: true,
                block_buffer: [0; BLOCK_SIZE],
            };
 

            // write inode to table here
            self.write_inode_to_table(0, &directory_inode);

            self.write_inode_to_table(file_table_entry.inode_num, &file_table_entry.inode);

            self.write_directory(&mut directory_contents, new_inode_num, filename);


        } else if let Some(existing_inode_num) = self.find_file_inode(filename) {

            // Use the current inode
            let mut existing_inode = self.get_inode_from_table(existing_inode_num);

            // Reset its stored-bytes and current-block to a state as if it were newly created
            existing_inode.bytes_stored = 0;

            // Clear the in-use bits for its existing data blocks, except for the first data block
            for i in 1..existing_inode.blocks.len() {
                let block_num = existing_inode.blocks[i] as usize;
                if block_num != 0 {
                    data_blocks_in_use[block_num / 8] &= !(1 << (block_num % 8));
                }
            }

            // Update the inode table with the modified inode
            self.write_inode_to_table(existing_inode_num, &existing_inode);

            // Create a file table entry for the existing file
            let file_table_entry = FileInfo {
                inode: existing_inode,
                inode_num: existing_inode_num,
                current_block: 0,
                offset: 0,
                writing: true,
                block_buffer: [0; BLOCK_SIZE],
            };

            // Add the file to the open table
            file_descriptor = self.add_to_open_table(file_table_entry).unwrap();


        } else {
            // Select an inode number and data block number
            let new_inode_num = Self::get_inode(&mut inode_block_buffer);
            let data_block_num = Self::get_inode(&mut data_blocks_in_use);
        
            // Create an inode for the new file
            let new_inode = Inode {
                bytes_stored: 0,
                blocks: [data_block_num as u8; MAX_FILE_BLOCKS],
            };
        
            // Save the new inode in the inode table
            self.write_inode_to_table(new_inode_num, &new_inode);

            self.write_directory(&mut directory_contents, new_inode_num, filename);
        
            // Create a file table entry for the newly created file
            let file_table_entry = FileInfo {
                inode: new_inode,
                inode_num: new_inode_num,
                current_block: 0,
                offset: 0,
                writing: true,
                block_buffer: [0; BLOCK_SIZE],
            };
        
            // Add the new file to the open table
            file_descriptor = self.add_to_open_table(file_table_entry).unwrap();

            self.write_inode_to_table(file_table_entry.inode_num, &file_table_entry.inode);
        
        }

        // Write back to the disk
        self.disk.write(0, &mut inode_block_buffer);
        self.disk.write(1, &mut data_blocks_in_use);
        
        // dont forget to write new_inode

        FileSystemResult::Ok(file_descriptor)

    }
    
    /* 
    pub fn write_itable(&mut self, inode_table: [u8; MAX_FILE_BYTES]) {
        let iterations = (inode_table.len() + self.block_buffer.len() - 1) / self.block_buffer.len();
        for i in 0..iterations {
            let start_index = i * self.block_buffer.len();
            let end_index = usize::min(inode_table.len(), start_index + self.block_buffer.len());
            let buffer_range = start_index..end_index;
    
            for (j, &inode_value) in inode_table[buffer_range].iter().enumerate() {
                if inode_value != 0 {
                    self.block_buffer[j] = inode_value;
                }
            }
            self.disk.write(i + 2, &mut self.block_buffer);
        }
    }
    */

    pub fn write_directory(&mut self, directory_contents: &mut [u8; MAX_FILE_BYTES], new_inode_num: usize, filename: &str) -> Result<(), ()> {
        let directory_inode: Inode<MAX_FILE_BLOCKS, BLOCK_SIZE> = self.get_inode_from_table(0);
    
        // Read directory contents from disk
        for i in 0..MAX_FILE_BLOCKS {
            if directory_inode.blocks[i] == 0 {
                break;
            }
            self.disk.read(directory_inode.blocks[i] as usize, &mut self.block_buffer);
            for j in 0..self.block_buffer.len() {
                directory_contents[i * self.block_buffer.len() + j] = self.block_buffer[j];
            }
        }
    
        // Write filename to directory_contents buffer
        for (i, &byte) in filename.as_bytes().iter().enumerate() {
            directory_contents[new_inode_num * MAX_FILENAME_BYTES + i] = byte;
        }
    
        // Write back to the disk
        for i in 0..MAX_FILE_BLOCKS {
            if directory_inode.blocks[i] == 0 {
                break;
            }
            for j in 0..self.block_buffer.len() {
                self.block_buffer[j] = directory_contents[i * self.block_buffer.len() + j];
            }
            self.disk.write(directory_inode.blocks[i] as usize, &mut self.block_buffer);
        }
    
        Ok(())
    }
    
    

    pub fn add_to_open_table(&mut self, file_info: FileInfo<MAX_FILE_BLOCKS, BLOCK_SIZE>) -> FileSystemResult<usize> {
        for (index, entry) in self.open.iter_mut().enumerate() {
            if entry.is_none() {
                *entry = Some(file_info);
                return FileSystemResult::Ok(index);
            }
        }
        FileSystemResult::Err(FileSystemError::TooManyOpen)
    }

    pub fn get_inode(block_buffer: &mut[u8; BLOCK_SIZE]) -> usize {
        for i in 0..block_buffer.len() {
            if block_buffer[i] < u8::MAX {
                let mut byte_value = block_buffer[i];
                let mut counter = 0;
                while byte_value % 2 == 1 {
                    byte_value >>= 1;
                    counter += 1
                }
                // set bit to 1
                block_buffer[i] |= 1 << counter; // correct?  
                return i * 8 + counter;

            }
        }
        // this is bad!!! 
        panic!("Disk full");
        // return FileSystemResult::Err(FileSystemError::DiskFull);
    }

    pub fn find_file_inode(&mut self, filename: &str) -> Option<usize> {
        let mut directory_contents = [0; MAX_FILE_BYTES];
        let mut temp_buffer = [0; BLOCK_SIZE];
    
        // Read the directory contents
        let directory_inode: Inode<MAX_FILE_BLOCKS, BLOCK_SIZE> = self.get_inode_from_table(0);

        for (i, &block) in directory_inode.blocks.iter().enumerate() {
            if block != 0 {
                self.disk.read(block as usize, &mut temp_buffer);
                directory_contents[i * BLOCK_SIZE..(i + 1) * BLOCK_SIZE].copy_from_slice(&temp_buffer);
            }
        }
    
        // Iterate through the directory contents and search for filename
        for i in 0..(MAX_FILE_BYTES / (MAX_FILENAME_BYTES + 2)) {
            let entry_start = i * (MAX_FILENAME_BYTES + 2);
            let entry_inode_num_start = entry_start + MAX_FILENAME_BYTES;
            let entry_inode_num = u16::from_le_bytes([directory_contents[entry_inode_num_start], directory_contents[entry_inode_num_start + 1]]) as usize;
    
            let mut entry_filename = String::new();
            for j in entry_start..entry_start + MAX_FILENAME_BYTES {
                if directory_contents[j] != 0 {
                    entry_filename.push(directory_contents[j] as char);
                } else {
                    break;
                }
            }
    
            if entry_filename == filename {
                return Some(entry_inode_num);
            }
        }
    
        None
    }
    
    

    pub fn get_inode_from_table<const MAX_BLOCKS: usize, const BLOCK_SIZ: usize>(
        &mut self,
        inode_num: usize,
    ) -> Inode<MAX_BLOCKS, BLOCK_SIZ> {
        // Calculate block number
        let inode_bytes = 2 + MAX_BLOCKS;
        let inodes_per_block = BLOCK_SIZ / inode_bytes;
        let block_num = INODE_TABLE_START + (inode_num / inodes_per_block);
        let offset = (inode_num % inodes_per_block) * inode_bytes;
    
        // Read the block containing inode
        self.disk.read(block_num, &mut self.block_buffer);
    
        // Extract inode information
        let bytes_stored = u16::from_le_bytes([ // I used this site but changed it to little endian https://stackoverflow.com/questions/29307474/how-can-i-convert-a-buffer-of-a-slice-of-bytes-u8-to-an-integer
            self.block_buffer[offset],
            self.block_buffer[offset + 1],
        ]);
        let mut blocks = [0; MAX_BLOCKS];
        blocks.copy_from_slice(&self.block_buffer[offset + 2..offset + 2 + MAX_BLOCKS]);
    
        Inode {
            bytes_stored,
            blocks,
        }
    }  

    pub fn write_inode_to_table(&mut self, inode_num: usize, inode: &Inode<MAX_FILE_BLOCKS, BLOCK_SIZE>) {
        // Calculate block and offset 
        let inode_block = INODE_TABLE_START + inode_num / self.inodes_per_block();
        let inode_offset = (inode_num % self.inodes_per_block()) * self.num_inode_bytes();


        // Read the block containing inode -- depends on if necessary
        self.disk.read(inode_block, &mut self.block_buffer);

        // Write the inode to the buffer
        self.block_buffer[inode_offset + 1] = (inode.bytes_stored >> 8) as u8;
        self.block_buffer[inode_offset] = (inode.bytes_stored &0xff) as u8; 
        // self.block_buffer[inode_offset + 2..inode_offset + 2 + MAX_FILE_BLOCKS] = inode.blocks; // may not work may have to do for loop
        for i in 0..MAX_FILE_BLOCKS {
            self.block_buffer[inode_offset + 2 + i] = inode.blocks[i];
        }        

        // Write updated block back to disk

        self.disk.write(inode_block, &mut self.block_buffer);
    }

    // I'm not using this and I probably should
    pub fn set_block_bit(buffer: &mut [u8], block_index: usize) {
        let byte_index = block_index / 8;
        let bit_index = block_index % 8;
        buffer[byte_index] |= 1 << bit_index;
    }

    fn allocate_data_block(&mut self) -> Result<[u8; 3], FileSystemError> {
        let mut data_table = [0; BLOCK_SIZE];
        self.disk.read(1, &mut data_table);
    
        let num_blocks = NUM_BLOCKS / 8;
        let mut block_allocated = false;
        let mut block_group = 0;
        let mut block_bit = 0;
        let mut block_num = 0;
    
        for (i, value) in data_table.iter().enumerate().take(num_blocks) {
            if *value != u8::MAX {
                for j in 0..8 {
                    if (value >> j) & 1 == 0 {
                        block_group = i;
                        block_bit = j;
                        block_num = i * 8 + j;
                        block_allocated = true;
                        break;
                    }
                }
            }
            if block_allocated {
                break;
            }
        }
    
        if !block_allocated {
            return Err(FileSystemError::DiskFull);
        }
    
        data_table[block_group] |= 1 << block_bit;
        self.disk.write(1, &data_table);
    
        Ok([block_group as u8, block_bit as u8, block_num as u8])
    }
    

    pub fn read(&mut self, fd: usize, buffer: &mut [u8]) -> FileSystemResult<usize> {
        if let Some(file_info) = self.open.get_mut(fd) {
            if let Some(file_info) = file_info.as_mut() { 
                let mut bytes_read = 0;
    
                while bytes_read < buffer.len() {
                    let current_block = file_info.current_block;
                    let current_offset = file_info.offset;
    
                    // Check if end of file has been reached
                    if bytes_read >= file_info.inode.bytes_stored as usize {
                        break;
                    }
    
                    // Read current block from disk into the block buffer
                    self.disk.read(file_info.inode.blocks[current_block] as usize, &mut self.block_buffer);
    
                    let bytes_left_in_block = BLOCK_SIZE - current_offset;
                    let bytes_to_copy = buffer.len() - bytes_read;
    
                    let bytes_to_copy = usize::min(bytes_left_in_block, bytes_to_copy);
                    let bytes_to_copy = usize::min(bytes_to_copy, file_info.inode.bytes_stored as usize - bytes_read);
    
                    // Copy bytes from the block buffer to user buffer
                    buffer[bytes_read..bytes_read + bytes_to_copy]
                        .copy_from_slice(&self.block_buffer[current_offset..current_offset + bytes_to_copy]);
    
                    bytes_read += bytes_to_copy;
    
                    // Check if the end of the current block has been reached
                    if current_offset + bytes_to_copy == BLOCK_SIZE {
                        file_info.current_block += 1; 
                        file_info.offset = 0; 
                    } else {
                        file_info.offset += bytes_to_copy;
                    }
                }
    
                FileSystemResult::Ok(bytes_read)
            } else {
                FileSystemResult::Err(FileSystemError::FileNotOpen)
            }
        } else {
            FileSystemResult::Err(FileSystemError::InvalidFileDescriptor(fd))
        }
    }
    
    pub fn write(&mut self, fd: usize, buffer: &[u8]) -> FileSystemResult<()> {
        if let Some(file_info) = self.open.get_mut(fd) {
            if let Some(file_info) = file_info {
                if !file_info.writing {
                    return FileSystemResult::Err(FileSystemError::NotOpenForWrite);
                }
                
                let mut bytes_written = 0;
                while bytes_written < buffer.len() {
                    let current_block = file_info.current_block;
                    let current_offset = file_info.offset;
        
                    // Copy bytes from the user's buffer to the block buffer
                    let bytes_to_copy = usize::min(BLOCK_SIZE - current_offset, buffer.len() - bytes_written);
    
                    file_info.block_buffer[current_offset..current_offset + bytes_to_copy]
                        .copy_from_slice(&buffer[bytes_written..bytes_written + bytes_to_copy]);
    
                    bytes_written += bytes_to_copy;
                    file_info.inode.bytes_stored += bytes_to_copy as u16;
    
                    // Check if the block buffer is full
                    if current_offset + bytes_to_copy == BLOCK_SIZE {
                        // Write the block buffer contents to the disk
                        self.disk.write(file_info.inode.blocks[current_block] as usize, &file_info.block_buffer);
    
                        file_info.current_block += 1;
                        file_info.offset = 0;
    
                        // Check if the file has reached max blocks
                        if file_info.current_block >= file_info.inode.blocks.len() {
                            return FileSystemResult::Err(FileSystemError::FileTooBig);
                        }
                    } else {
                        file_info.offset += bytes_to_copy;
                    }
                }
    
                FileSystemResult::Ok(())
            } else {
                FileSystemResult::Err(FileSystemError::FileNotOpen)
            }
        } else {
            FileSystemResult::Err(FileSystemError::InvalidFileDescriptor(fd))
        }
    }
    
    pub fn add_to_inode(&mut self, inode_num: u8, new_data_block : u8) -> [u8;MAX_FILE_BYTES] {
        let mut buffer = [0; BLOCK_SIZE];
        let inode_start = inode_num as u16 * self.num_inode_bytes() as u16;

        let mut last_block = 0;

        let mut index = (inode_start + 3) as u16;
        let mut flag = false;

        for i in (inode_start + 2) as u16..(inode_start as u16 + self.num_inode_bytes() as u16){
            if !flag {

                flag = true;
                last_block = self.file_content_buffer[i as usize];

            } else if flag && self.file_content_buffer[i as usize] == last_block{

                break;

            } else if flag {

                last_block = self.file_content_buffer[i as usize];
                index = i as u16;

            }
        }
        self.file_content_buffer[index as usize] = new_data_block;
        self.write_itable(2 + inode_start as usize / BLOCK_SIZE)
        
    }
    
    pub fn write_itable(&mut self, start_block: usize)  -> [u8;MAX_FILE_BYTES]{
        let mut count = 0;
        let mut inode_buffer = [0; BLOCK_SIZE];
        let mut start = start_block;

        let total = MAX_FILES_STORED * self.num_inode_bytes();
        let mut total_count = 0;

        for i in self.file_content_buffer { 

            if total_count >= total {
                break;
            }

            if count + 1 == BLOCK_SIZE {

                self.disk.write(start, &mut inode_buffer);

                inode_buffer = [0; BLOCK_SIZE];
                count = 0;
                start += 1;

            } else {

                inode_buffer[count] = i;

            }
            
            count += 1;
            total_count += 1;
        }
        return self.file_content_buffer;
    }
    
    pub fn close(&mut self, fd: usize) -> FileSystemResult<()> {
        if let Some(file_info_option) = self.open.get_mut(fd) {
            if let Some(mut file_info) = file_info_option.take() { // Remove the entry in the file table
                if file_info.writing {
                    // Update the inode to store its new size
                    let new_size = (file_info.current_block * BLOCK_SIZE) + file_info.offset;
                    file_info.inode.bytes_stored = new_size as u16;
    
                    // Write the updated inode back to the disk
                    let inode_block = INODE_TABLE_START + file_info.inode_num / self.inodes_per_block();
                    let inode_offset = (file_info.inode_num % self.inodes_per_block()) * self.num_inode_bytes();
                    self.disk.read(inode_block, &mut self.block_buffer);
                    self.block_buffer[inode_offset] = (new_size >> 8) as u8;
                    self.block_buffer[inode_offset + 1] = new_size as u8;
                    self.disk.write(inode_block, &self.block_buffer);
                }
    
                // Mark the inode as closed in open_inodes
                self.open_inodes[file_info.inode_num] = false;
                FileSystemResult::Ok(())

            } else {
                FileSystemResult::Err(FileSystemError::FileNotOpen)
            }
        } else {
            FileSystemResult::Err(FileSystemError::InvalidFileDescriptor(fd))
        }
    }
    

}


//              :'(



#[cfg(test)]
mod tests {
    use super::*;

    const BLOCK_SIZE: usize = 64;
    const MAX_FILES_STORED: usize = 32;

    fn make_small_fs() -> FileSystem<16, 64, 255, 8, 512, 32, 8> {
        FileSystem::new(ramdisk::RamDisk::new())
    }

    #[test]
    fn test_short_write() {
        let mut sys = make_small_fs();
        let f1 = sys.open_create("one.txt").unwrap();
        sys.assert_block(0, 0, &[3, 0]);
        sys.assert_block(1, 0, &[255, 1, 0]);
        sys.assert_block(2, 0, &[16, 0, 7]);
        sys.assert_block(2, 10, &[0, 0, 8]);
        sys.assert_block(7, 0, &[0, 0, 0, 0, 0, 0, 0, 0, 111, 110, 101, 46, 116, 120, 116, 0]);
        sys.write(f1, "This is a test.".as_bytes()).unwrap();
        let mut buffer = [0; 50];
        sys.close(f1).unwrap();
        sys.assert_block(8, 0, &[84, 104, 105, 115, 32, 105, 115, 32, 97, 32, 116, 101, 115, 116, 46]);
        sys.assert_block(2, 0, &[16, 0, 7]);
        sys.assert_block(2, 10, &[15, 0, 8]);
        let f2 = sys.open_read("one.txt").unwrap();
        let bytes_read = sys.read(f2, &mut buffer).unwrap();
        assert_eq!(bytes_read, 15);
        let s = core::str::from_utf8(&buffer[0..bytes_read]).unwrap();
        assert_eq!(s, "This is a test.");
    }
        
    const LONG_DATA: &str = "This is a much, much longer message.
    It crosses a number of different lines in the text editor, all synthesized
    with the goal of exceeding the 64 byte block limit by a considerable amount.
    To that end, this text contains considerable excessive verbiage.";

    #[test]
    fn test_long_write() {
        assert_eq!(265, LONG_DATA.len());
        let mut sys = make_small_fs();
        let f1 = sys.open_create("one.txt").unwrap();
        sys.write(f1, LONG_DATA.as_bytes()).unwrap();
        sys.close(f1);
        sys.assert_block(0, 0, &[3, 0, 0]);
        sys.assert_block(1, 0, &[255, 31, 0]);
        sys.assert_block(2, 0, &[16, 0, 7]);
        sys.assert_block(2, 10, &[9, 1, 8, 9, 10, 11, 12]);
        let read = read_to_string(&mut sys, "one.txt");
        assert_eq!(read.as_str(), LONG_DATA);
    }

    fn read_to_string(
        sys: &mut FileSystem<16, BLOCK_SIZE, 255, 8, 512, 32, 8>,
        filename: &str,
    ) -> String {
        let fd = sys.open_read(filename).unwrap();
        let mut read = String::new();
        let mut buffer = [0; 10];
        loop {
            let num_bytes = sys.read(fd, &mut buffer).unwrap();
            let s = core::str::from_utf8(&buffer[0..num_bytes]).unwrap();
            read.push_str(s);
            if num_bytes < buffer.len() {
                sys.close(fd).unwrap();
                return read;
            }
        }
    }

    #[test]
    fn test_complex_1() {
        let one = "This is a message, a short message, but an increasingly long message.
        This is a message, a short message, but an increasingly long message.";
        let two = "This is the second message I have chosen to undertake in this particular test.
        This is a continuation of this ever-so-controversial second message.\n";
        let mut sys = make_small_fs();
        let f1 = sys.open_create("one.txt").unwrap();
        sys.write(f1, one[0..one.len() / 2].as_bytes()).unwrap();
        let f2 = sys.open_create("two.txt").unwrap();
        sys.write(f2, two[0..two.len() / 2].as_bytes()).unwrap();
        sys.write(f1, one[one.len() / 2..one.len()].as_bytes())
            .unwrap();
        sys.write(f2, two[two.len() / 2..two.len()].as_bytes())
            .unwrap();
        sys.close(f1).unwrap();
        sys.close(f2).unwrap();
        assert_eq!(one, read_to_string(&mut sys, "one.txt").as_str());
        assert_eq!(two, read_to_string(&mut sys, "two.txt").as_str());
    }

    // #[test]
    // fn test_complex_2() {
    //     let one = "This is a message, a short message, but an increasingly long message.
    //     This is a message, a short message, but an increasingly long message.";
    //     let two = "This is the second message I have chosen to undertake in this particular test.
    //     This is a continuation of this ever-so-controversial second message.\n";
    //     let mut sys = make_small_fs();
    //     let f1 = sys.open_create("one.txt").unwrap();
    //     sys.write(f1, one[0..one.len() / 2].as_bytes()).unwrap();
    //     let f2 = sys.open_create("two.txt").unwrap();
    //     sys.write(f2, two[0..two.len() / 2].as_bytes()).unwrap();
    //     sys.close(f1).unwrap();
    //     sys.close(f2).unwrap();

        // let f3 = sys.open_append("two.txt").unwrap();
        // let f4 = sys.open_append("one.txt").unwrap();
    //     sys.write(f4, one[one.len() / 2..one.len()].as_bytes())
    //         .unwrap();
    //     sys.write(f3, two[two.len() / 2..two.len()].as_bytes())
    //         .unwrap();
    //     sys.close(f1).unwrap();
    //     sys.close(f2).unwrap();
    //     assert_eq!(one, read_to_string(&mut sys, "one.txt").as_str());
    //     assert_eq!(two, read_to_string(&mut sys, "two.txt").as_str());
    // }
    
    #[test]
    fn test_complex_3() {
        let one = "This is a message, a short message, but an increasingly long message.
        This is a message, a short message, but an increasingly long message.";
        let two = "This is the second message I have chosen to undertake in this particular test.
        This is a continuation of this ever-so-controversial second message.\n";
        let mut sys = make_small_fs();
        let f1 = sys.open_create("one.txt").unwrap();
        sys.write(f1, one.as_bytes()).unwrap();
        sys.close(f1).unwrap();

        let f2 = sys.open_create("one.txt").unwrap();
        sys.write(f2, two.as_bytes()).unwrap();
        sys.close(f2).unwrap();

        assert_eq!(two, read_to_string(&mut sys, "one.txt").as_str());
    }

    #[test]
    fn test_file_not_found() {
        let mut sys = make_small_fs();
        let f1 = sys.open_create("one.txt").unwrap();
        sys.write(f1, "This is a test.".as_bytes()).unwrap();
        sys.close(f1).unwrap();
        match sys.open_read("one.tx") {
            FileSystemResult::Ok(_) => panic!("Shouldn't have found the file"),
            FileSystemResult::Err(e) => assert_eq!(e, FileSystemError::FileNotFound),
        }
    }

    #[test]
    fn test_file_not_open() {
        let mut sys = make_small_fs();
        let f1 = sys.open_create("one.txt").unwrap();
        sys.write(f1, "This is a test.".as_bytes()).unwrap();
        sys.close(f1).unwrap();
        let fd = sys.open_read("one.txt").unwrap();
        let mut buffer = [0; 10];
        match sys.read(fd + 1, &mut buffer) {
            FileSystemResult::Ok(_) => panic!("Should be an error!"),
            FileSystemResult::Err(e) => assert_eq!(e, FileSystemError::FileNotOpen),
        }
    }

    #[test]
    fn test_not_open_for_read() {
        let mut sys = make_small_fs();
        let f1 = sys.open_create("one.txt").unwrap();
        sys.write(f1, "This is a test.".as_bytes()).unwrap();
        let mut buffer = [0; 10];
        match sys.read(f1, &mut buffer) {
            FileSystemResult::Ok(_) => panic!("Should not work!"),
            FileSystemResult::Err(e) => assert_eq!(e, FileSystemError::NotOpenForRead),
        }
    }

    #[test]
    fn test_not_open_for_write() {
        let mut sys = make_small_fs();
        let f1 = sys.open_create("one.txt").unwrap();
        sys.write(f1, "This is a test.".as_bytes()).unwrap();
        sys.close(f1).unwrap();
        let f2 = sys.open_read("one.txt").unwrap();
        match sys.write(f2, "this is also a test".as_bytes()) {
            FileSystemResult::Ok(_) => panic!("Should be an error"),
            FileSystemResult::Err(e) => assert_eq!(e, FileSystemError::NotOpenForWrite),
        }
    }

    #[test]
    fn test_filename_too_long() {
        let mut sys = make_small_fs();
        match sys.open_create("this_is_an_exceedingly_long_filename_to_use.txt") {
            FileSystemResult::Ok(_) => panic!("This should be an error"),
            FileSystemResult::Err(e) => assert_eq!(e, FileSystemError::FilenameTooLong),
        }
    }

    #[test]
    fn test_already_open() {
        let mut sys = make_small_fs();
        let f1 = sys.open_create("one.txt").unwrap();
        sys.write(f1, "This is a test.".as_bytes()).unwrap();
        match sys.open_read("one.txt") {
            FileSystemResult::Ok(_) => panic!("Should be an error"),
            FileSystemResult::Err(e) => assert_eq!(e, FileSystemError::AlreadyOpen),
        }
    }

    #[test]
    fn test_file_too_big() {
        let mut sys = make_small_fs();
        let f1 = sys.open_create("one.txt").unwrap();
        for _ in 0..sys.max_file_size() - 1 {
            sys.write(f1, "A".as_bytes()).unwrap();
        }
        match sys.write(f1, "B".as_bytes()) {
            FileSystemResult::Ok(_) => panic!("Should be an error!"),
            FileSystemResult::Err(e) => assert_eq!(e, FileSystemError::FileTooBig),
        }
    }

    #[test]
    fn test_too_many_files() {
        let mut sys = make_small_fs();
        for i in 0..MAX_FILES_STORED - 1 {
            let filename = format!("file{i}");
            let f = sys.open_create(filename.as_str()).unwrap();
            let content = format!("This is sentence {i}");
            sys.write(f, content.as_bytes()).unwrap();
            sys.close(f).unwrap();
        }
        match sys.open_create("Final") {
            FileSystemResult::Ok(_) => panic!("This should be an error!"),
            FileSystemResult::Err(e) => assert_eq!(e, FileSystemError::TooManyFiles),
        }
    }
}