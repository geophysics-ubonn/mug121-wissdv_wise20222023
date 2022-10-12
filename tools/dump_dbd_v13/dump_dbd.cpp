/*
dump_dbd - read dataTaker binary data files

Copyright (c) 2005-2010 Thermo Fisher Scientific Australia Pty. Ltd.

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <ctype.h>

#ifdef __GNUC__
#include <limits.h>
#else
#include <io.h>
#endif

#define VERSION "13"

#ifdef __GNUC__
#if (LONG_MAX > 0x7fffffff)
typedef unsigned long   uint64;     // Linux 64-bit, GNU g++
#else
typedef unsigned long long uint64;  // Linux 32-bit, GNU g++
#endif
typedef unsigned int    uint32;
typedef int             int32;
typedef unsigned short  uint16;
typedef unsigned char   uint8;
#else
typedef unsigned __int64 uint64;    // WIN32, Visual C++
typedef unsigned long   uint32;
typedef long            int32;
typedef unsigned short  uint16;
typedef unsigned char   uint8;
#endif

#define TRUE        true
#define FALSE       false

// REV32 and REV16 macros assume that dump_dbd is running in a little-endian (Intel) environment

#define TAG(a,b,c,d) ( ((a)<<24) | ((b)<<16) | ((c)<<8) | (d) )
#define REV32(a)     ( (a) = ( (((a)&0x000000ff)<<24) | (((a)&0x0000ff00)<<8) | (((a)&0x00ff0000)>>8) | (((a)&0xff000000)>>24) ) )
#define REV16(a)     ( (a) = ( (((a)&0x00ff)<<8) | (((a)&0xff00)>>8) ) )

#define DTDT_TAG    TAG('D','T','D','T')
#define HDR1_TAG    TAG('H','D','R','1')
#define DHDR_TAG    TAG('D','H','D','R')
#define AHDR_TAG    TAG('A','H','D','R')
#define DATA_TAG    TAG('D','A','T','A')
#define ALRM_TAG    TAG('A','L','R','M')
#define EOF_TAG     TAG('E','O','F','!')

enum
{
    PREFIX_DATATYPE_FLOATING_POINT  = 0x10,
    PREFIX_DATATYPE_INTEGER         = 0x20,
    PREFIX_DATATYPE_ASCII_TEXT      = 0x30,
    PREFIX_DATATYPE_ABSOLUTE_TIME   = 0x40,
    PREFIX_DATATYPE_RELATIVE_TIME   = 0x50,
    PREFIX_DATATYPE_DATA_RECORD     = 0x60,
    PREFIX_DATATYPE_ALARM_RECORD    = 0x70,
    PREFIX_DATATYPE_COMPACT_DATA    = 0x80,
};

enum
{
    PREFIX_DATASTATE_VALID_DATA         = 0x01,
    PREFIX_DATASTATE_NOT_YET_SET        = 0x02,
    PREFIX_DATASTATE_DATA_OVER_RANGE    = 0x03,
    PREFIX_DATASTATE_DATA_UNDER_RANGE   = 0x04,
    PREFIX_DATASTATE_OPEN_CIRCUIT       = 0x06,
    PREFIX_DATASTATE_REFERENCE_ERROR    = 0x07,
    PREFIX_DATASTATE_CANNOT_STORE       = 0x08,
    PREFIX_DATASTATE_CANNOT_REPRESENT   = 0x09,
    PREFIX_DATASTATE_DISCONTINUITY      = 0x0b,
};

enum
{
    PREFIX_LENGTH_MASK              = 0xff00,
    PREFIX_DATATYPE_MASK            = 0x00f0,
    PREFIX_DATASTATE_MASK           = 0x000f
};

#define PREFIX_LENGTH(a) (((uint32)(a) & PREFIX_LENGTH_MASK) >> 8)
#define PREFIX_TYPE(a)    ((uint32)(a) & PREFIX_DATATYPE_MASK)
#define PREFIX_STATE(a)   ((uint32)(a) & PREFIX_DATASTATE_MASK)

#pragma pack(1)

typedef struct
{
    uint32 tag;
    uint32 length;
}
SectionPrefix;

typedef struct
{
    char            userChanName[16];
    char            defaultChanName[16];
    char            units[16];
    char            unused1[1];
    char            displayFormat[7];
    char            userChanNameExt[8];
    uint16          chanNumber;
    uint16          prefix;
    char            unused2[12];
}
ChannelMapEntry;



typedef union
{
    float           fp;
    uint32          i32;
    char            text[1];
    uint8           time[8];
}
DataValue;

typedef struct
{
    uint16          prefix;
    uint32          timestampOffset;
    DataValue       data;
}
DataPoint;

typedef struct
{
    uint16          prefix;
    DataValue       data;
}
CompactDataPoint;

typedef struct
{
    uint16          prefix;
    uint8           tstamp[8];
    DataPoint       dataPoint[1];
}
DataRecord;

typedef struct
{
    uint16          prefix;
    uint8           tstamp[8];
    uint8           alarmNumber;
    uint8           alarmState;
    char            alarmText[245];
}
AlarmRecord;

typedef struct
{
    uint32          bytesPerRecord;
    uint32          oldestRecord;
    uint32          newestRecord;
}
RecordInfo;

//----file sections

// HDR1
typedef struct
{
    char            jobName[16];
    char            sha1[20];
    uint32          schedule;
    int32           timezone;
    char            schedName[16];
    char            serialNum[12];
}
HeaderSection;

// DHDR
typedef struct
{
    RecordInfo      recInfo;
    uint32          dataPointsPerRecord;
    ChannelMapEntry channelMap[1];  // dataPointsPerRecord elements
}
DataHeaderSection;

// AHDR
typedef struct
{
    RecordInfo      recInfo;
}
AlarmHeaderSection;

#pragma pack()


typedef struct
{
    HeaderSection*      header;
    DataHeaderSection*  dataHeader;
    AlarmHeaderSection* alarmHeader;
    DataRecord*         dataRec;       // storage for one record
    AlarmRecord*        alarmRec;
} Storefile;


static bool readDataRec( uint32 rec );
static bool readAlarmRec( uint32 rec );
static bool readSectionPrefix( SectionPrefix* pprefix );
static bool readBlock( void* buf, uint32 size );
static void formatTimestamp( uint64 timestamp, const char* fmt, char* buf );
static void formatInterval( uint64 timestamp, const char* fmt, char* buf );

static char **makeFileList( char **argv, int first_file_arg_index, int argc );

static FILE* fp;
static Storefile file;
static bool verbose;
static bool showData;
static bool showAlarms;
static bool showHeaders;
static bool showColNames;
static bool useTimeOffsets;
static bool ignoreHeaderFormat;
static char columnNames[32768];

static void rev64( uint8* p )
{
    uint8 t;
    t=p[0]; p[0]=p[7]; p[7]=t;
    t=p[1]; p[1]=p[6]; p[6]=t;
    t=p[2]; p[2]=p[5]; p[5]=t;
    t=p[3]; p[3]=p[4]; p[4]=t;
}

int main(int argc, char* argv[])
{
    int fileArg = argc;
    long pos, eof;
    int dtdtNum;
    SectionPrefix secPrefix;
    RecordInfo* recInfo;
    bool isData;
    int i;
    char* pn;
    int rv=2;
    char *path = NULL;
    int fileNum = 1;
    char  **filesToDump = NULL;
    int   list_index = 0;

    for( i=1; i<argc; i++ )
    {
        if( argv[i][0]=='-' )
        {
            switch( argv[i][1] )
            {
            case 'v': verbose = TRUE; break;
            case 'd': showData = TRUE; break;
            case 'a': showAlarms = TRUE; break;
            case 'h': showHeaders = TRUE; break;
            case 'n': showColNames = TRUE; break;
            case 'i': ignoreHeaderFormat = TRUE; break;
            case 't': useTimeOffsets = TRUE; break; // not implemented
            }
        }
        else
        {
            fileArg = i;
            break;
        }
    }
    if( fileArg == argc )
    {
        fprintf( stderr, "\ndump_dbd Ver " VERSION " - dump contents of specified .DBD file(s)\n"
                         " Usage:    dump_dbd [options] <infile>\n"
                         " Options:  -v  verbose\n"
                         "           -h  show headers only\n"
                         "           -d  show data only\n"
                         "           -a  show alarms only\n"
                         "           -n  show column names (headings)\n"
                         "           -i  ignore format spec in header\n"

                         );
        return 1;
    }
    if( !showData && !showAlarms && !showHeaders )
    {
        showData = showAlarms = showHeaders = TRUE;
    }

    /*
     * Make a list of DBD files to dump.  The list is derived from the command
     * line arguments.
     */
    filesToDump = makeFileList( argv, fileArg, argc );

    for ( list_index = 0; filesToDump[ list_index ]; ++list_index ) 
    {
        bool firstDtDt = TRUE;
        path = filesToDump[ list_index ];
        eof = -1;
        pos = 0;
        dtdtNum = 0;

        if( (fp=fopen( path, "rb" )) == 0 )
        {
            printf( "\n? Cannot open %s for read\n", path );
            goto nextfile;
        }

        if( showHeaders )
            printf( "\nFile name: %s\n", path );
        else
            printf( "\n" );

        for(;;)
        {
            // read section prefix
            bool gotsect = readSectionPrefix( &secPrefix );
            pos = ftell( fp );
            if( firstDtDt && (!gotsect || (secPrefix.tag != DTDT_TAG)) )
            {
                printf( "? Cannot read DTDT section\n" );
                goto nextfile;
            }
            firstDtDt = FALSE;
            if( !gotsect )
            {
                // reached physical eof
                goto nextfile;
            }
            if( (eof >= 0) && (pos > eof) && (secPrefix.tag != DTDT_TAG) )
            {
                // past eof specified in previous DTDT, and there is not another DTDT
                goto nextfile;
            }
            pos = ftell( fp );

            // read remainder of section based on section type
            switch( secPrefix.tag )
            {
            case DTDT_TAG:
                dtdtNum++;
                if( showHeaders ) printf( "\nDTDT #%d start: 0x%lx size: 0x%lx\n", dtdtNum, pos, secPrefix.length );
                eof = pos + secPrefix.length;
                break;

            case HDR1_TAG:
                if( showHeaders ) printf( "\nHDR1 start: 0x%lx size: 0x%lx\n", pos, secPrefix.length );
                if( (secPrefix.length < sizeof(HeaderSection)) || (secPrefix.length > 0x100000) )
                {
                    printf("?%08x Invalid HDR1 length\n",secPrefix.length);
                    fseek( fp, secPrefix.length, SEEK_CUR );
                    break;
                }
                file.header = (HeaderSection *)malloc( secPrefix.length );
                if( ! readBlock( file.header, secPrefix.length ) )
                {
                    goto nextfile;
                }
                REV32( file.header->schedule );
                REV32( file.header->timezone );
                if( showHeaders ) 
                {
                    printf( "Job:         %.16s\n", file.header->jobName );
                    printf( "Schedule:    %c  \"%.16s\"\n", file.header->schedule, file.header->schedName );
                    printf( "Serial num:  %.12s\n", file.header->serialNum );
                    printf( "Timezone:    " );
                    if( file.header->timezone == -1 )
                    {
                        printf( "none\n" );
                    }
                    else if( (file.header->timezone > -24*60) && (file.header->timezone < 24*60) )
                    {
                        printf( "%+d:%02d\n", file.header->timezone/60, file.header->timezone%60 );
                    }
                    else
                    {
                        printf( "?%08x Bad timezone\n", file.header->timezone );
                    }
                    printf( "SHA-1:       " );
                    for( i=0; i<20; i++ )
                    {
                        printf( "%02x", (uint8)file.header->sha1[i] );
                    }
                    printf( "\n" );
                }
                break;

            case DHDR_TAG:
            case AHDR_TAG:
                isData = (secPrefix.tag == DHDR_TAG);
                if( showHeaders ) printf( "\n%s start: 0x%lx size: 0x%lx\n", isData ? "DHDR" : "AHDR", pos, secPrefix.length );
                if( (secPrefix.length < (isData ? sizeof(DataHeaderSection) : sizeof(AlarmHeaderSection))) || (secPrefix.length > 0x100000) )
                {
                    printf("?%08x Invalid %s length\n", secPrefix.length, isData ? "DHDR" : "AHDR");
                    fseek( fp, secPrefix.length, SEEK_CUR );
                    break;
                }
                if( isData )
                {
                    file.dataHeader = (DataHeaderSection *)malloc( secPrefix.length );
                    if( ! readBlock( file.dataHeader, secPrefix.length ) ) goto nextfile;
                    recInfo = &file.dataHeader->recInfo;
                }
                else
                {
                    file.alarmHeader = (AlarmHeaderSection *)malloc( secPrefix.length );
                    if( ! readBlock( file.alarmHeader, secPrefix.length ) ) goto nextfile;
                    recInfo = &file.alarmHeader->recInfo;
                }

                REV32( recInfo->bytesPerRecord );
                REV32( recInfo->oldestRecord );
                REV32( recInfo->newestRecord );
                if( showHeaders ) 
                {
                    printf( "Oldest record:      0x%x\n", recInfo->oldestRecord );
                    printf( "Newest record:      0x%x\n", recInfo->newestRecord );
                    printf( "Bytes/record:       0x%x\n", recInfo->bytesPerRecord );
                }

                if( isData )
                {
                    REV32( file.dataHeader->dataPointsPerRecord );
                    if( showHeaders ) printf( "DataPoints/record:  0x%x\n", file.dataHeader->dataPointsPerRecord );
                    if( secPrefix.length < sizeof(DataHeaderSection)+(file.dataHeader->dataPointsPerRecord-1)*sizeof(ChannelMapEntry) )
                    {
                        printf( "? DHDR length not consistent with datapts/rec value\n" );
                        break;
                    }
                    if( showHeaders ) printf( "Channel map:\n" );
                    pn = columnNames;
                    pn += sprintf( pn, "\"Timestamp\"" );
                    uint32 recbytes = sizeof(DataRecord)-sizeof(DataPoint);
                    for( i=0; i < (int)file.dataHeader->dataPointsPerRecord; i++ )
                    {
                        ChannelMapEntry* cmap = &file.dataHeader->channelMap[ i ];

                        REV16( cmap->chanNumber );
                        REV16( cmap->prefix );
                        if( showHeaders ) printf( "  %2u: %.16s  name:\"%.16s%.8s\"  units:\"%.16s\"  fmt:\"%.8s\"  prefix:0x%04hx\n",
                            i, cmap->defaultChanName, cmap->userChanName, cmap->userChanNameExt, cmap->units, 
                            (cmap->unused1[0] ? cmap->unused1 : cmap->displayFormat), cmap->prefix );

                        pn += sprintf( pn, ",\"%.16s%.8s", cmap->userChanName, cmap->userChanNameExt );
                        if( *cmap->units ) pn += sprintf( pn, " (%.16s)", cmap->units );
                        pn += sprintf( pn, "\"" );
                        recbytes += 2 + PREFIX_LENGTH(cmap->prefix);
                    }
                    pn += sprintf( pn, "\n" );
                    if( recbytes > recInfo->bytesPerRecord )
                    {
                        printf( "? CMAP length not consistent with bytes/rec value\n" );
                        break;
                    }
                }
                break;

            case DATA_TAG:
            case ALRM_TAG:
                isData = (secPrefix.tag == DATA_TAG);
                if( showHeaders ) printf( "\n%s start: 0x%lx size: 0x%lx\n", isData ? "DATA" : "ALRM", pos, secPrefix.length );
                recInfo = isData ? &file.dataHeader->recInfo : &file.alarmHeader->recInfo;
                if( (isData && !showData) || (!isData && !showAlarms) )
                {
                    // skip to next section
                    fseek( fp, secPrefix.length, SEEK_CUR );
                }
                else if( (int)recInfo->oldestRecord < 0 )
                {
                    // no records: skip to next section
                    if( verbose && ((isData && showData) || (!isData && showAlarms)) ) printf( "No records\n" );
                    fseek( fp, secPrefix.length, SEEK_CUR );
                }
                else
                {
                    uint32 oldest, newest;
                    uint32 logged, capacity;
                    uint32 dataStartPos;
                    uint32 rec, num;

                    if( showColNames ) printf( isData ? columnNames : "\"Timestamp\",\"Num\",\"Transition\",\"Text\"\n" );

                    // allocate a buffer to hold one record
                    if( isData )
                        file.dataRec = (DataRecord *)malloc( recInfo->bytesPerRecord );
                    else
                        file.alarmRec = (AlarmRecord *)malloc( recInfo->bytesPerRecord );
                    capacity = secPrefix.length / recInfo->bytesPerRecord;
                    oldest = recInfo->oldestRecord;
                    newest = recInfo->newestRecord;
                    logged = (oldest <= newest) ? (newest - oldest + 1) : (capacity - (oldest - newest - 1));
                    if( verbose && ((isData && showData) || (!isData && showAlarms)) ) printf( "%u records\n", logged );
                    // seek to oldest record
                    dataStartPos = ftell(fp);
                    fseek( fp, dataStartPos + (oldest * recInfo->bytesPerRecord), SEEK_SET );
                    // read data records from oldest to newest
                    for( rec = oldest, num = 0; num < logged; num++ )
                    {
                        // read and display one record
                        if( isData )
                            readDataRec( rec );
                        else
                            readAlarmRec( rec );

                        // next record, wrap if required
                        if( ++rec == capacity )
                        {
                            rec = 0;
                            fseek( fp, dataStartPos, SEEK_SET );
                        }
                    }
                    // seek to end of section
                    fseek( fp, dataStartPos + secPrefix.length, SEEK_SET );
                }
                break;

            case EOF_TAG:
                if( verbose && showHeaders ) printf( "\nEOF! start: 0x%lx size: 0x%lx\n", pos, secPrefix.length );
                break;

            default:
                printf( "\n?%04x Unknown section tag (at 0x%lx)\n", secPrefix.tag, pos );
                goto nextfile;
            }
        }
nextfile:
        free( file.header ); file.header = 0;
        free( file.dataHeader ); file.dataHeader = 0;
        free( file.alarmHeader );  file.alarmHeader = 0;
        free( file.dataRec ); file.dataRec = 0;
        free( file.alarmRec ); file.alarmRec = 0;

        if ( fp ) fclose(fp);
#ifndef __GNUC__
        free( path );
#endif
    }

    rv = 0;
    if( verbose ) printf( "--End--\n" );

    free( filesToDump );
    return rv;
}


static bool readDataRec( uint32 rec )
{
    uint64 timestamp;
    uint32 i;
    DataPoint* dp;
    char buf[60];
    bool isCompact=FALSE;
    bool isValid=FALSE;
    long pos = ftell(fp);

    if( ! readBlock( file.dataRec, file.dataHeader->recInfo.bytesPerRecord ) )
    {
        return FALSE;
    }

    // output timestamp (also rec number and record prefix if verbose)
    // abort record output if invalid prefix type is not "data record"
    REV16( file.dataRec->prefix );
    switch( PREFIX_TYPE(file.dataRec->prefix) )
    {
    case PREFIX_DATATYPE_COMPACT_DATA: isCompact = isValid = TRUE; break;
    case PREFIX_DATATYPE_DATA_RECORD: isValid = TRUE; break;
    }
    if( !isValid )
    {
        if( verbose ) printf( "%5u: ", rec );
        printf( "?%04hx Bad record prefix (at 0x%x)\n", file.dataRec->prefix, pos );
        return TRUE;
    }
    if( verbose ) printf( "%5u: %04hx,", rec, file.dataRec->prefix );
    rev64( file.dataRec->tstamp );
    timestamp = *(uint64 *)file.dataRec->tstamp;
    formatTimestamp( timestamp, "datetime", buf );
    printf( buf );

    // output datapoints
    for( i=0, dp = &file.dataRec->dataPoint[0]; i < file.dataHeader->dataPointsPerRecord; i++ )
    {
        const char* dispFmt = ignoreHeaderFormat ? "" : file.dataHeader->channelMap[ i ].displayFormat;
        uint32 dplen;
        char* p;
        int txtlen;
        DataValue* pdv;

        REV16( dp->prefix );
        // only the length part of the data record prefix needs to match that in channel map
        if( (dp->prefix & PREFIX_LENGTH_MASK) == (file.dataHeader->channelMap[ i ].prefix & PREFIX_LENGTH_MASK) )
        {
            pdv = isCompact ? &((CompactDataPoint *)dp)->data : &dp->data;
            REV32( dp->timestampOffset );
            dplen = PREFIX_LENGTH(dp->prefix);
            printf( ", " );
            if( verbose ) 
            {
                printf( "%04hx,", dp->prefix );
                if( !isCompact ) 
                    printf( "T+%u.%04u,", dp->timestampOffset/1000000, ((dp->timestampOffset%1000000) + 50) / 100 );
            }

            switch( PREFIX_STATE(dp->prefix) )
            {
            case PREFIX_DATASTATE_VALID_DATA:
                switch( PREFIX_TYPE(dp->prefix) )
                {
                case PREFIX_DATATYPE_FLOATING_POINT:
                    {
                    char fmt[5] = "%.7g";
                    REV32( pdv->i32 );
                    // standard display format [FEM]n
                    if( dispFmt[0] && strchr( "FEM", dispFmt[0] ) && isdigit( dispFmt[1] ) )
                    {
                        fmt[2] = dispFmt[1];
                        switch( dispFmt[0] )
                        {
                        case 'F': fmt[3] = 'f'; break;
                        case 'E': fmt[3] = 'E'; break;
                        case 'M': fmt[3] = 'G'; break;
                        }
                    }
                    // legacy display format %.n[fEG] (first char not part of dispFmt any more)
                    else if( (dispFmt[0] == '.') && isdigit( dispFmt[1] ) && strchr( "fEG", dispFmt[2] ) )
                    {
                        fmt[2] = dispFmt[1];
                        fmt[3] = dispFmt[2];
                    }
                    printf( fmt, pdv->fp );
                    }
                    break;
                case PREFIX_DATATYPE_INTEGER:
                    REV32( pdv->i32 );
                    printf( "%d", pdv->i32 );
                    break;
                case PREFIX_DATATYPE_ASCII_TEXT:
                    // text length = datapoint length length less 4 bytes for timestamp offset if present
                    txtlen = PREFIX_LENGTH(file.dataRec->prefix) - (isCompact ? 4 : 0);
                    // replace control chars with spaces
                    for( p=pdv->text; *p && (p < &pdv->text[txtlen]); p++ )
                    {
                        if( *p < ' ' ) *p = ' ';
                    }
                    printf( "\"%.*s\"", txtlen, pdv->text );
                    break;
                case PREFIX_DATATYPE_ABSOLUTE_TIME:
                    rev64( pdv->time );
                    timestamp = *(uint64 *)pdv->time;
                    formatTimestamp( timestamp, dispFmt, buf );
                    printf( buf );
                    break;
                case PREFIX_DATATYPE_RELATIVE_TIME:
                    rev64( pdv->time );
                    timestamp = *(uint64 *)pdv->time;
                    formatInterval( timestamp, 0, buf );
                    printf( buf );
                    break;
                default:
                    printf( "?%08x", pdv->i32 );
                    break;
                }
                break;

            case PREFIX_DATASTATE_NOT_YET_SET:
                printf( "NotYetSet" );
                break;

            case PREFIX_DATASTATE_DATA_OVER_RANGE:
                printf( "Overrange" );
                break;

            case PREFIX_DATASTATE_DATA_UNDER_RANGE:
                printf( "Underrange" );
                break;

            case PREFIX_DATASTATE_REFERENCE_ERROR:
                printf( "RefError" );
                break;

            case PREFIX_DATASTATE_DISCONTINUITY:
                printf( "Discontinuity" );
                break;

            default:
                printf( "?%08x", pdv->i32 );
                break;
            }
        }
        else
        {
            // prefix doesn't match that in cmap: length may be invalid so don't read any more datapoints in this record
            printf( "?%04x Bad datapoint prefix", dp->prefix );
            break;
        }

        // advance to next datapoint
        dp = (DataPoint*) ((uint8*)dp + (((dp->prefix & PREFIX_LENGTH_MASK) >> 8) + 2));
    }
    printf( "\n" );
    return TRUE;
}

static bool readAlarmRec( uint32 rec )
{
    uint64 timestamp;
    char* p;
    char buf[60];
    int txtlen;
    long pos = ftell(fp);

    if( ! readBlock( file.alarmRec, file.alarmHeader->recInfo.bytesPerRecord ) )
    {
        return FALSE;
    }

    // output timestamp (also rec number and record prefix if verbose)
    // abort record output if invalid prefix type is not "alarm record"
    REV16( file.alarmRec->prefix );
    if( (PREFIX_TYPE(file.alarmRec->prefix) != PREFIX_DATATYPE_ALARM_RECORD) ||
        ((2 + PREFIX_LENGTH(file.alarmRec->prefix)) != file.alarmHeader->recInfo.bytesPerRecord) )
    {
        if( verbose ) printf( "%5u: ", rec );
        printf( "?%04hx Bad record prefix (at 0x%x)\n", file.alarmRec->prefix, pos );
        return TRUE;
    }
    if( verbose ) printf( "%5u: %04hx,", rec, file.alarmRec->prefix );
    rev64( file.alarmRec->tstamp );
    timestamp = *(uint64 *)file.alarmRec->tstamp;
    formatTimestamp( timestamp, "datetime", buf );
    printf( buf );

    // text length = record length length less 10 bytes for timestamp, alarm number, alarm transition
    txtlen = PREFIX_LENGTH(file.alarmRec->prefix) - 8 - 1 - 1;
    // replace control chars with spaces
    for( p=file.alarmRec->alarmText; *p && (p < &file.alarmRec->alarmText[txtlen]); p++ )
    {
        if( *p < ' ' ) *p = ' ';
    }
    printf( ",ALARM%u,%u,\"%.*s\"\n", 
        file.alarmRec->alarmNumber, file.alarmRec->alarmState, txtlen, file.alarmRec->alarmText );

    return TRUE;
}

static bool readSectionPrefix( SectionPrefix* pprefix )
{
    if( fread( pprefix, 1, sizeof(SectionPrefix), fp ) == sizeof(SectionPrefix) )
    {
        REV32( pprefix->tag );
        REV32( pprefix->length );
        return TRUE;
    }
    else if( ferror(fp) )
    {
        perror( "? Error reading file" );
        return FALSE;
    }
    else
    {
        return FALSE;
    }
}

static bool readBlock( void* buf, uint32 size )
{
    long pos = ftell(fp);
    uint32 numread = fread( buf, 1, size, fp );

    if( numread == size )
    {
        return TRUE;
    }
    else if( ferror(fp) )
    {
        perror( "? Error reading file" );
        return FALSE;
    }
    else
    {
        printf( "? Cannot read file (at 0x%x, read 0x%x bytes, expected 0x%x)\n", pos, numread, size );
        return FALSE;
    }
}

static void formatTimestamp( uint64 timestamp, const char* fmt, char* buf )
{
    time_t sec = (time_t)(timestamp/1000000); // sec since 1/1/1989
    uint32 subsec = ((uint32)(timestamp%1000000) + 50) / 100;
    struct tm* tm;


    sec += ((1989-1970)*365 + 5) * 24*60*60;  // sec since 1/1/1970
    tm = gmtime( &sec );
    if( tm )
    {
        // standard display format: T or D
        // legacy display format: time or date or datetime (first char not part of disp_fmt any more)
        if( (fmt[0] == 'T') || (fmt[0] == 'i') )
        {
            sprintf( buf, "%02d:%02d:%02d.%04d",
                tm->tm_hour, tm->tm_min, tm->tm_sec, subsec );
        }
        else if( (fmt[0] == 'D') || ((fmt[0] == 'a') && (fmt[3] == 0)) )
        {
            sprintf( buf, "%u/%02d/%02d",
                tm->tm_year+1900, tm->tm_mon+1, tm->tm_mday );
        }
        else
        {
            sprintf( buf, "%u/%02d/%02d %02d:%02d:%02d.%04d",
                tm->tm_year+1900, tm->tm_mon+1, tm->tm_mday, tm->tm_hour, tm->tm_min, tm->tm_sec, subsec );
        }
    }
    else
    {
        sprintf( buf, "?TS:%08x-%08x", timestamp>>32, timestamp&0xffffffff );
    }
}

static void formatInterval( uint64 timestamp, const char* fmt, char* buf )
{
    time_t sec = (time_t)(timestamp/1000000); // sec
    uint32 subsec = ((uint32)(timestamp%1000000) + 50) / 100;
    struct tm* tm;

    tm = gmtime( &sec );
    sprintf( buf, "%d:%02d:%02d.%04d",
        tm->tm_yday*24 + tm->tm_hour, tm->tm_min, tm->tm_sec, subsec );
}


// Return an array of ptrs to pathnames
// for Linux, these ptrs point to argv strings
// for Windows, these ptrs point to dynamically allocated memory, which needs to be freed
// In both cases the array itself is dynamically allocated memory, which needs to be freed
static char **makeFileList( char **argv, int file_arg_index, int argc )
{
    enum { LIST_SZ_INCREMENT = 0x10 };

    char **list = NULL;
    int  list_index = 0;
    int  list_sz = LIST_SZ_INCREMENT;

#ifdef __GNUC__
    // for Linux, shell will already have expanded any command line wildcards
    list = ( char ** )calloc( 1 + ( argc - file_arg_index ), sizeof( char * ) );

    for ( ; file_arg_index < argc; ++file_arg_index, ++list_index )
    {
        list[ list_index ] = argv[ file_arg_index ];
    }

    list[ list_index ] = NULL;
    
#else
    // for Windows, we need to process wildcards ourselves
    struct _finddata_t  find_info;
    long                find_handle;

    char drive[_MAX_DRIVE];
    char dir[_MAX_DIR];
    char fname[_MAX_FNAME];
    char ext[_MAX_EXT];
 
    list = ( char ** )calloc( LIST_SZ_INCREMENT, sizeof( char * ) );

    while ( file_arg_index < argc )
    {
        find_handle = _findfirst( argv[ file_arg_index ], &find_info );
        if( find_handle == -1L )
        {
            fprintf( stderr, "dump_dbd: no files matching '%s'\n", argv[file_arg_index] );
            goto done;
        }

        do
        {
            if ( find_info.attrib & _A_SUBDIR ) continue;
            _splitpath( argv[ file_arg_index ], drive, dir, fname, ext );

            if ( list_index >= ( list_sz - 1 ) )
            {
                list_sz += LIST_SZ_INCREMENT;
                list = ( char ** )realloc( list, list_sz * sizeof( char * ) );
            }

            list[ list_index ] = (char *)malloc( _MAX_PATH );

            sprintf( list[ list_index ], "%s%s%s", drive, dir, find_info.name );

            list[ ++list_index ] = NULL;
        }
        while( _findnext( find_handle, &find_info ) == 0 );
        ++file_arg_index;
    }

done:
    if( find_handle != -1L ) _findclose( find_handle );

#endif

    return list;
}

