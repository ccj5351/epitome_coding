#include "makeDirectory.h"

// to create a directory
// returns true if successful;
// otherwise, returns false;
bool MakeDir(const string & dir_path){
	bf::path p(dir_path); 
	if (bf::exists(p) && is_directory(p)){ // if directory p actually exist;
		std::cout << p << " already exists.\n";
		return true;
	}
	else // p actually does not exist;
		return bf::create_directory(p);
}


bool find_file( const bf::path & dir_path,     // in this directory,
                const std::string & file_name, // search for this name,
                bf::path & path_found          // placing path here if found
){
  if ( !bf::exists( dir_path ) ) return false;
  bf::directory_iterator end_itr; // default construction yields past-the-end
  for ( bf::directory_iterator itr( dir_path );
        itr != end_itr;
        ++itr ){
    if ( bf::is_directory(itr->status())){
      if ( find_file( itr->path(), file_name, path_found ) ) return true;
    }
	else if ( itr->path().filename() == file_name ){ // see below
      path_found = itr->path();
      return true;
    }
  }
  return false;
}


/*
// see http://www.cplusplus.com/forum/general/61834/
As follows: - from boost path.hpp
#   ifdef BOOST_WINDOWS_API
    const std::string string() const { return string(codecvt()); } 
    const std::string string(const codecvt_type& cvt) const
    { 
      std::string tmp;
      if (!m_pathname.empty())
        path_traits::convert(&*m_pathname.begin(), &*m_pathname.begin()+m_pathname.size(),
          tmp, cvt);
      return tmp;
    }
    
    //  string_type is std::wstring, so there is no conversion
    const std::wstring&  wstring() const { return m_pathname; }
    const std::wstring&  wstring(const codecvt_type&) const { return m_pathname; }
	*/


void GetDirList(const string & directory, vector<string> * dirlist){
	 bf::path p(directory);
  if (bf::is_directory(p)){
    for (bf::directory_iterator itr(p); itr!=bf::directory_iterator(); ++itr){
		if (bf::is_directory(itr->status()))
			dirlist -> push_back(itr->path().filename().string());
    }
  }
  else 
	  cout <<  p << " is not a directory!\n";
}


// to get the names of all the image files in the current path
void GetFileList(const string& directory, vector<string>* filelist){
	bf::path p(directory);
  if (bf::is_directory(p)){
    for (bf::directory_iterator itr(p); itr!=bf::directory_iterator(); ++itr){
		if (bf::is_regular_file(itr->status()))
			filelist -> push_back(itr->path().filename().string());
    }
  }
  else 
	  cout <<  p << " is not a directory!\n";
}

// to get the names of all the image files in the current path
void GetFileList(const string& directory, 
	const string & fileExtension, // e.g., = ".txt", including the dot ".";
	vector<string>* filelist
){
	bf::path p(directory);
  if (bf::is_directory(p)){
    for (bf::directory_iterator itr(p); itr!=bf::directory_iterator(); ++itr){
		if ((bf::is_regular_file(itr->status())) && (itr->path().extension().string() == fileExtension))
			filelist -> push_back(itr->path().filename().string());
    }
  }
  else 
	  cout <<  p << " is not a directory!\n";
}