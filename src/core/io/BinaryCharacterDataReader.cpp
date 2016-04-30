#include "BinaryCharacterDataReader.h"

#include "ContinuousBinaryTaxonData.h"
#include "DiscreteBinaryTaxonData.h"
#include "Exception.h"
#include "FileManager.h"
#include "StringUtilities.h"

/** Returns whether a file exists */
bool BinaryCharacterDataReader::fileExists(const char* fn) const {
    
	bool exists = false;
	FILE *fp = fopen(fn, "r");
	if (fp != NULL)
    {
		fclose(fp);
		exists = true;
    }
	return exists;
}


std::string BinaryCharacterDataReader::findFileNameFromPath(const std::string& fp) const {
    
    std::string::size_type pos = fp.find_last_of('/');
    if ( pos != std::string::npos )
    {
        std::string fn = fp.substr(pos+1, fp.size()-pos-1);
        return fn;
    }
    return "";
}


/** Get a reference to this singleton object */
BinaryCharacterDataReader& BinaryCharacterDataReader::getInstance(void) {
    
	static BinaryCharacterDataReader rb;
	return rb;
}



/** Format the error exception string for problems specifying the file/path name */
void BinaryCharacterDataReader::formatError(FileManager& fm, std::string& errorStr) {
    
    bool fileNameProvided    = fm.isFileNamePresent();
    bool isFileNameGood      = fm.testFile();
    bool isDirectoryNameGood = fm.testDirectory();
    
    if ( fileNameProvided == false && isDirectoryNameGood == false )
        errorStr += "Could not read contents of directory \"" + fm.getFilePath() + "\" because the directory does not exist";
    else if (fileNameProvided == true && (isFileNameGood == false || isDirectoryNameGood == false))
    {
        errorStr += "Could not read file named \"" + fm.getFileName() + "\" in directory named \"" + fm.getFilePath() + "\" ";
        if (isFileNameGood == false && isDirectoryNameGood == true)
            errorStr += "because the file does not exist";
        else if (isFileNameGood == true && isDirectoryNameGood == false)
            errorStr += "because the directory does not exist";
        else
            errorStr += "because neither the directory nor the file exist";
    }
}


/** Try to determine if the file is likely to be in Fasta format */
bool BinaryCharacterDataReader::isFastaFile(std::string& fn) {
    
    // open file
    std::ifstream fStrm;
    fStrm.open(fn.c_str(), std::ios::in);
    
    // read the file token-by-token looking for Fasta things
    int ch = fStrm.get();
    fStrm.unget();
    std::string line;
    int lineNum = 0, lastCarotLine = -100;
    int numSequences = 0;
    
    while (getline(fStrm,line))
    {
        std::string word;
        std::istringstream lStrm(line);
        
        int wordNum = 0;
        
        bool comment = false;
        
        while(lStrm >> word)
        {
            if(comment)
                continue;
            
            if (wordNum == 0 && word[0] == '>')
            {
                if (lineNum - lastCarotLine > 1)
                {
                    lastCarotLine = lineNum;
                    numSequences++;
                }
                else{
                    std::cerr << lineNum << "\t" << lastCarotLine << std::endl;
                    return false;
                }
            }
            else if (wordNum == 0 && word[0] == ';')
            {
                comment = true;
            }
            else if (lineNum > lastCarotLine && word[0] != '>' && word[0] != ';')
            {
                for(size_t i=0; i<word.size(); i++)
                {
                    if(word[i] != '0' && word[i] != '1' && word[i] != '?' && word[i] != '-')
                        return false;
                }
            }
            
            wordNum++;
        }
        
        lineNum++;
    }
    
    // close file
    fStrm.close();
    
    if (numSequences < 1)
    {
        std::cerr << numSequences << " < 1" << std::endl;
        return false;
    }
    
    return true;
}

/** Try to determine if the file is likely to be in Fasta format */
bool BinaryCharacterDataReader::isPastaFile(std::string& fn) {
    
    // open file
    std::ifstream fStrm;
    fStrm.open(fn.c_str(), std::ios::in);
    
    // read the file token-by-token looking for Fasta things
    int ch = fStrm.get();
    fStrm.unget();
    std::string line;
    int lineNum = 0, lastCarotLine = -100;
    int numSequences = 0;
    while (getline(fStrm,line))
    {
        std::string word;
        std::istringstream lStrm(line);
        
        int wordNum = 0;
        bool comment = false;
        while(lStrm >> word)
        {
            if(comment)
                continue;
            
            if (wordNum == 0 && word[0] == '>')
            {
                if (lineNum - lastCarotLine > 1)
                {
                    lastCarotLine = lineNum;
                    numSequences++;
                }
                else{
                    std::cerr << lineNum << "\t" << lastCarotLine << std::endl;
                    return false;
                }
            }
            else if (wordNum == 0 && word[0] == ';')
            {
                comment = true;
            }
            else if (lineNum > lastCarotLine && word[0] != '>' && word[0] != ';')
            {
                if(IsFloat(word))
                {
                    RealNumber f = atof(word.c_str());
                    
                    if(f > 1.0 || f < 0.0)
                    {
                        //std::cerr << word << " notin (0,1)" << std::endl;
                        return false;
                    }
                }
                else if(word != "?" && word != "-"){
                    //std::cerr << "'" << word << "' not float or gap" << std::endl;
                    return false;
                }
            }
            
            wordNum++;
        }
        
        lineNum++;
    }
    
    // close file
    fStrm.close();
    
    if (numSequences < 1)
    {
        std::cerr << numSequences << " < 1" << std::endl;
        return false;
    }
    
    return true;
}


/**
 * Try to determine if the file is likely to be in Nexus format. We check if the first word in the
 * file is #NEXUS. If not, then we check if the file name ending is ".nex". If neither is true, it
 * is probably not a NEXUS file.
 */
bool BinaryCharacterDataReader::isNexusFile(const std::string& fn) {
    
    // open file, read first word, close file
	std::ifstream fStrm;
    fStrm.open(fn.c_str(), std::ios::in);
    std::string word;
    fStrm >> word;
    fStrm.close();
    
    if (word=="#NEXUS")
        return true;
    else {
        size_t found = fn.find_last_of(".");
        if ( found != std::string::npos && fn.substr(found+1) == "nex" )
            return true;
        else
            return false;
    }
}


/** Try to determine if the file is likely to be in Phylip format */
bool BinaryCharacterDataReader::isPhylipFile(std::string& fn, bool& isInterleaved) {
    
    // open file
	std::ifstream fStrm;
    fStrm.open(fn.c_str(), std::ios::in);
    std::string seqStr = "";
    
    // read the file token-by-token looking for NEXUS things
    bool foundNumTaxa = false, foundNumChar = false;
    unsigned int numTaxa = 0;
    std::vector<std::string> taxonNames;
    int ch = fStrm.get();
    fStrm.unget();
    std::string word = "";
    int wordNum = 0, lineNum = 0;
    while (ch != EOF)
    {
        word = "";
        fStrm >> word;
        StringUtilities::toLower( word );
        
        if (lineNum == 0 && wordNum == 0 && StringUtilities::isNumber(word) == true)
        {
            std::istringstream buf(word);
            buf >> numTaxa;
            foundNumTaxa = true;
        }
        else if (lineNum == 0 && wordNum == 1 && StringUtilities::isNumber(word) == true)
            foundNumChar = true;
        else if (lineNum > 0 && wordNum == 0 && word != "" && word.size() < 100)
            taxonNames.push_back( word );
        else if (lineNum > 0 && wordNum > 0)
            seqStr += word;
        
        wordNum++;
        ch = fStrm.get();
        if (ch == '\n' || ch == '\r' || ch == EOF)
        {
            lineNum++;
            wordNum = 0;
        }
    }
    
    // close file
    fStrm.close();
    
    if (foundNumTaxa == true && foundNumChar == true)
    {
        if (taxonNames.size() == 0)
            return false;
        if (taxonNames.size() % numTaxa != 0)
            return false;
        
        if (taxonNames.size() > numTaxa)
            isInterleaved = true;

        return true;
    }
    
    return false;
}


BinaryCharacterData* BinaryCharacterDataReader::readMatrix(std::string fn) {
    
    // check that the file/path name has been correctly specified
    FileManager myFileManager( fn );
    if ( myFileManager.getFileName() == "" && myFileManager.getFilePath() == "" )
    {
        std::string errorStr = "";
        formatError(myFileManager, errorStr);
        throw Exception("Could not find file or path with name \"" + fn + "\"");
    }
    
    // set up a vector of strings containing the name or names of the files to be read
    
    // Set up a map with the file name to be read as the key and the file type as the value. Note that we may not
    // read all of the files in the string called "vectorOfFileNames" because some of them may not be in a format
    // that can be read.
    BinaryCharacterData* m = NULL;
    bool isInterleaved = false;
    bool pasta = false;
    if (isNexusFile(fn) == true)
        m = ReadNexus(fn);
    else if (isPhylipFile(fn, isInterleaved))
    {
        if (isInterleaved == true)
            m = ReadPhylip(fn,1);
        else
            m = ReadPhylipSequential(fn);
    }
    else if (isPastaFile(fn))
    {
        pasta = true;
        m = ReadPasta(fn);
    }
    else if (isFastaFile(fn))
    {
        m = ReadFasta(fn);
    }
    else
        addWarning("Unknown file type");
    
    // print summary of results of file reading to the user
    if (m != NULL)
    {
        std::cout << "Successfully read file" << std::endl;
        if(pasta)
            std::cout << "Read " << ((ContinuousBinaryCharacterData*)m)->getNumberOfCharacters() << " characters, " << m->getNumberOfTaxa() << " taxa" << std::endl;
        else
            std::cout << "Read " << m->getNumberOfCharacters() << " characters, " << m->getNumberOfTaxa() << " taxa" << std::endl;
    }
    else
    {
        std::set<std::string> myWarnings = getWarnings();
        if ( myWarnings.size() > 0 )
        {
            std::stringstream o3;
            o3 << "Did not read the file for the following ";
            if (myWarnings.size() == 1)
                o3 << "reason:";
            else
                o3 << "reasons:";
            std::cerr << o3.str() << std::endl;
            for (std::set<std::string>::iterator it = myWarnings.begin(); it != myWarnings.end(); it++)
                std::cerr << "* "+(*it) << std::endl;
        }
    }
    
    return m;
}

BinaryCharacterData* BinaryCharacterDataReader::ReadNexus(std::string filespec)   {

    std::ifstream theStream(filespec.c_str());
    
    std::vector<DiscreteBinaryTaxonData*> Data;
    
    try {

        GoPastNextWord(theStream, "dimensions");
        GoPastNext(theStream, '=');
        size_t Ntaxa;
        theStream >> Ntaxa;
        GoPastNext(theStream, '=');
        size_t Nsite;
        theStream >> Nsite;
        GoPastNextWord(theStream, "format");
        GoPastNextWord(theStream, "datatype");
        GoPastNext(theStream, '=');
        std::string type;
        theStream >> type;
        


        if (!EquivalentStrings(type,"restriction") && !EquivalentStrings(type,"standard") && !EquivalentStrings(type,"binary"))  {
            std::cerr << "error cannot recognize data type\n";
            std::cerr << type << "\n";
            exit(1);
        }
        
        for(size_t i =0; i < Ntaxa; i++)
        {
            DiscreteBinaryTaxonData* td = new DiscreteBinaryTaxonData;
            Data.push_back(td);
        }

        GoPastNextWord(theStream, "Matrix");

        size_t l = 0;
        while (l<Nsite) {
            size_t m = 0;
            for (size_t i=0; i<Ntaxa; i++) {
                std::string temp;
                theStream >> temp;
                while (temp == "[") {
                    unsigned char c;
                    c = 'i';
                    while (c != ']') c = theStream.get();
                    theStream >> temp;
                }

                if (!l) {
                    Data[i]->setTaxonName(temp);
                }
                else    {
                    if (temp != Data[i]->getTaxonName())    {
                        std::cerr << "error when reading tree base: " << temp << '\t' << Data[i]->getTaxonName() << '\n';
                        exit(1);
                    }
                }

                unsigned char c;
                size_t k = l;
                do  {
                    c = theStream.get();
                    if (c == '[')   {
                        while (c != ']') c = theStream.get();
                        c = theStream.get();
                    }
                    if ((c != ' ') && (c != '\t') && (c != '\n') && (c != 13))  {
                        Data[i]->addCharacter(false);
                        if (c == '(')   {
                            Data[i]->setGap(k,true);
                            while (c != ')')    {
                                theStream >> c;
                            }
                        }
                        else if (c == '{')  {
                            Data[i]->setGap(k,true);
                            while (c != '}')    {
                                theStream >> c;
                            }
                        }
                        else if (c == '-' || c == '?')  {
                            Data[i]->setGap(k,true);
                        }
                        else    {
                            std::ostringstream s;
                            s << c;
                            bool val;
                            std::istringstream(s.str()) >> val;
                            Data[i]->setCharacter(k,val);
                        }
                        k++;
                    }
                }
                while ((!theStream.eof()) && (c != '\n') && (c != 13));
                if (theStream.eof())    {
                    if (i < Ntaxa-1)    {
                        std::cerr << "error : found " << i << " taxa instead of " << Ntaxa << " in datafile\n";
                        exit(1);
                    }
                }
                if (!m) {
                    m = k;
                }
                else    {
                    if (m != k) {
                        std::cerr << "error when reading nexus : " << m << '\t' << k << '\n';
                        std::cerr << "taxa : " << i << '\t' << Data[i]->getTaxonName() << '\n';
                        if (m > k)  {
                            while (k != m)  {
                                Data[i]->setGap(k,true);
                                k++;
                            }
                        }
                    }
                }
            }
            l= m;
        }
    }
    catch(...)  {
        std::cerr << "error while reading data file\n";
        return 0;
    }
    
    BinaryCharacterData rt;
    
    for(size_t i = 0; i < Data.size(); i++)
        rt.addTaxonData(Data[i]);
    
    return rt.clone();
}

BinaryCharacterData* BinaryCharacterDataReader::ReadPhylipSequential(std::string filespec)  {

    std::ifstream theStream(filespec.c_str());
    
    BinaryCharacterData Data;
    
    try {

        // std::cerr << "beware: phylip data sets only for amino acids for the moment\n";
        // std::cerr.flush();
    
        std::string temp;
        theStream >> temp;
        if (!IsInt(temp))   {
            std::cerr << "error when reading data\n";
            std::cerr << "data should be formatted as follows:\n";
            std::cerr << "#taxa #sites\n";
            std::cerr << "name1 seq1.....\n";
            std::cerr << "name2 seq2.....\n";
            std::cerr << "...\n";
            std::cerr << '\n';
            exit(1);
        }
        size_t Ntaxa = Int(temp);
        theStream >> temp;
        if (!IsInt(temp))   {
            std::cerr << "error when reading data\n";
            std::cerr << "data should be formatted as follows:\n";
            std::cerr << "#taxa #sites\n";
            std::cerr << "name1 seq1.....\n";
            std::cerr << "name2 seq2.....\n";
            std::cerr << "...\n";
            std::cerr << '\n';
            exit(1);
        }
        size_t Nsite = Int(temp);

        size_t ntaxa = 0;
        while ((!theStream.eof()) && (ntaxa<Ntaxa)) {
            DiscreteBinaryTaxonData* td = new DiscreteBinaryTaxonData;
            
            theStream >> temp;
            td->setTaxonName(temp);
            size_t nsite = 0;

            char c;
            do  {
                c = theStream.get();
                if ((!theStream.eof()) && (c != ' ') && (c != '\n') && (c != '\t') && (c!=13))  {
                    td->addCharacter(false);
                    if (c == '(')   {
                        td->setGap(nsite,true);
                        while (c != ')')    {
                            theStream >> c;
                        }
                    }
                    else if (c == '{')  {
                        td->setGap(nsite,true);
                        while (c != '}')    {
                            theStream >> c;
                        }
                    }
                    else if (c == '-' || c == '?')  {
                        td->setGap(nsite,true);
                    }
                    else    {
                        std::ostringstream s;
                        s << c;
                        bool val;
                        std::istringstream(s.str()) >> val;
                        td->setCharacter(nsite,val);
                    }
                    nsite++;
                }
            }
            while ((!theStream.eof()) && (nsite < Nsite));
            ntaxa++;
            
            Data.addTaxonData(td);
        }
    }
    catch(...)  {
        std::cerr << "error while reading data file\n";
        exit(1);
    }
    
    return Data.clone();
}

BinaryCharacterData* BinaryCharacterDataReader::ReadFasta(std::string filespec)  {

    std::ifstream theStream(filespec.c_str());
        
    BinaryCharacterData Data;
    
    try {

        size_t ntaxa = 0;
        std::string line;
        
        size_t nsite = 0;
        int lineNum = 0, lastCarotLine = 0;
        
        DiscreteBinaryTaxonData* td;
        
        while (getline(theStream,line)) {
            std::string word;
            std::istringstream lStrm(line);
            
            int wordNum = 0;           
            bool comment = false;
 
            while(lStrm >> word)
            {
                if(comment)
                    continue;
                
                if (wordNum == 0 && word[0] == '>')
                {
                    if(ntaxa > 0)
                        Data.addTaxonData(td);
                    
                    lastCarotLine = lineNum;
                    ntaxa++;
                    
                    word.erase(0, 1);
                    
                    td = new DiscreteBinaryTaxonData(word);
                    
                    nsite = 0;
                }
                else if (word[0] == ';')
                {
                    comment = true;
                }
                else if (lineNum > lastCarotLine)
                {
                    std::istringstream wStrm(word);
                    char c;
                    while(wStrm.get(c)){
                        td->addCharacter(false);
                        if (c == '(')   {
                            td->setGap(nsite,true);
                            while (c != ')')    {
                                wStrm >> c;
                            }
                        }
                        else if (c == '{')  {
                            td->setGap(nsite,true);
                            while (c != '}')    {
                                wStrm >> c;
                            }
                        }
                        else if (c == '-' || c == '?')  {
                            td->setGap(nsite,true);
                        }
                        else    {
                            std::ostringstream s;
                            s << c;
                            bool val;
                            std::istringstream(s.str()) >> val;
                            td->setCharacter(nsite,val);
                        }
                        nsite++;
                    }
                }
                
                wordNum++;
            }
            
            lineNum++;
        }
        
        Data.addTaxonData(td);
    }
    catch(...)  {
        std::cerr << "error while reading data file\n";
        exit(1);
    }
                
    return Data.clone();
}

BinaryCharacterData* BinaryCharacterDataReader::ReadPhylip (std::string filespec, int repeattaxa) {

    std::ifstream theStream(filespec.c_str());
    
    std::vector<BinaryTaxonData *> Data;
    
    try {

        std::string temp;
        theStream >> temp;
        if (!IsInt(temp))   {
            std::cerr << "error when reading data\n";
            std::cerr << "data should be formatted as follows:\n";
            std::cerr << "#taxa #sites\n";
            std::cerr << "name1 seq1.....\n";
            std::cerr << "name2 seq2.....\n";
            std::cerr << "...\n";
            std::cerr << '\n';
            exit(1);
        }
        size_t Ntaxa = Int(temp);
        theStream >> temp;
        if (!IsInt(temp))   {
            std::cerr << "error when reading data\n";
            std::cerr << "data should be formatted as follows:\n";
            std::cerr << "#taxa #sites\n";
            std::cerr << "name1 seq1.....\n";
            std::cerr << "name2 seq2.....\n";
            std::cerr << "...\n";
            std::cerr << '\n';
            exit(1);
        }
        size_t Nsite = Int(temp);
        std::cerr << Ntaxa << '\t' << Nsite << '\n';

        for(size_t i = 0; i < Ntaxa; i++)
        {
            DiscreteBinaryTaxonData* td = new DiscreteBinaryTaxonData;
            Data.push_back(td);
        }
        int l = 0;
        int block = 0;
        while (l<Nsite) {
            block++;
            int m = 0;
            for (int i=0; i<Ntaxa; i++) {
                if ((!l) || repeattaxa) {
                    std::string temp;
                    theStream >> temp;
                    if (!l) {
                        Data[i]->setTaxonName(temp);
                    }
                    else    {
                        if (temp != Data[i]->getTaxonName())    {
                            std::cerr << "error when reading data: read " << temp << " instead of " << Data[i]->getTaxonName() << '\n';
                            exit(1);
                        }
                    }
                }

                unsigned char c;
                int k = l;
                do  {

                    c = theStream.get();
                    if ((!theStream.eof()) && (c != ' ') && (c != '\n') && (c != '\t') && (c!=13))  {
                        if (c == '(')   {
                            Data[i]->setGap(k,true);
                            while (c != ')')    {
                                theStream >> c;
                            }
                        }
                        else if (c == '{')  {
                            Data[i]->setGap(k,true);
                            while (c != '}')    {
                                theStream >> c;
                            }
                        }
                        else if (c == '-' || c == '?')  {
                            Data[i]->setGap(k,true);
                        }
                        else    {
                            std::ostringstream s;
                            s << c;
                            bool val;
                            std::istringstream(s.str()) >> val;
                            Data[i]->setCharacter(k,val);
                        }
                        k++;
                    }
                }
                while ((!theStream.eof()) && (c != '\n') && (c != 13));
                if (theStream.eof())    {
                    if (i < Ntaxa-1)    {
                        std::cerr << "error : found " << i << " taxa instead of " << Ntaxa << " in datafile\n";
                        exit(1);
                    }
                }
                c = theStream.peek();
                while ((!theStream.eof()) && ((c == '\n') || (c == 13)))    {
                    c = theStream.get();
                    c = theStream.peek();
                }
                
                if (!m) {
                    m = k;
                }
                else    {
                    if (m != k) {
                        std::cerr << "error when reading data non matching number of sequences in block number " << block << " for taxon " << i << " " << Data[i]->getTaxonName() << '\n';
                        std::cerr << "taxa : " << i << '\t' << Data[i]->getTaxonName() << '\n';
                        std::cerr << "read " << k << " instead of " << m << "characters\n";
                        exit(1);
                    }
                }
            }
            l= m;
        }
        if (l<Nsite)    {
            std::cerr << "error : reached end of stream \n";
            std::cerr << "data should be formatted as follows:\n";
            std::cerr << "#taxa #sites\n";
            std::cerr << "name1 seq1.....\n";
            std::cerr << "name2 seq2.....\n";
            std::cerr << "...\n";
            std::cerr << '\n';
            exit(1);
        }
    }
    catch(...)  {
        std::cerr << "error while reading data file\n";
    }
    
    BinaryCharacterData rt;
    
    for(size_t i = 0; i < Data.size(); i++)
        rt.addTaxonData(Data[i]);
    
    return rt.clone();
}

ContinuousBinaryCharacterData* BinaryCharacterDataReader::ReadPasta(std::string filespec)  {

    std::ifstream theStream(filespec.c_str());
    
    ContinuousBinaryCharacterData Data;
    
    try {

        size_t ntaxa = 0;
        std::string line;
        
        size_t nsite = 0;
        int lineNum = 0, lastCarotLine = 0;
        
        ContinuousBinaryTaxonData* td;
        
        while (getline(theStream,line)) {
            std::string word;
            std::istringstream lStrm(line);
            
            int wordNum = 0;           
            bool comment = false;
 
            while(lStrm >> word)
            {
                if(comment)
                    continue;
                
                if (wordNum == 0 && word[0] == '>')
                {
                    if(ntaxa > 0)
                        Data.addTaxonData(td);
                    
                    lastCarotLine = lineNum;
                    ntaxa++;
                    
                    word.erase(0, 1);
                    
                    td = new ContinuousBinaryTaxonData(word);
                    
                    nsite = 0;
                }
                else if (word[0] == ';')
                {
                    comment = true;
                }
                else if (lineNum > lastCarotLine)
                {
                    td->addCharacter(1.0);
                    if(word == "?" || word == "-")
                        td->setGap(nsite,true);
                    else
                        td->setCharacter(nsite, atof(word.c_str()));
                    
                    nsite++;
                }
                
                wordNum++;
            }
            
            lineNum++;
        }
        
        Data.addTaxonData(td);
    }
    catch(...)  {
        std::cerr << "error while reading data file\n";
        exit(1);
    }
                
    return Data.clone();
}
