/* =========================================================================================

   This is an auto-generated file: Any edits you make may be overwritten!

*/

#pragma once

namespace BinaryData
{
    extern const char*   lmao_png;
    const int            lmao_pngSize = 5829447;

    extern const char*   lilypad_png;
    const int            lilypad_pngSize = 5866916;

    extern const char*   pinkpad_png;
    const int            pinkpad_pngSize = 5982915;

    extern const char*   Catalogue_2_0_ttf;
    const int            Catalogue_2_0_ttfSize = 27948;

    extern const char*   stems_png;
    const int            stems_pngSize = 4968865;

    extern const char*   catfader_png;
    const int            catfader_pngSize = 3970235;

    // Number of elements in the namedResourceList and originalFileNames arrays.
    const int namedResourceListSize = 6;

    // Points to the start of a list of resource names.
    extern const char* namedResourceList[];

    // Points to the start of a list of resource filenames.
    extern const char* originalFilenames[];

    // If you provide the name of one of the binary resource variables above, this function will
    // return the corresponding data and its size (or a null pointer if the name isn't found).
    const char* getNamedResource (const char* resourceNameUTF8, int& dataSizeInBytes);

    // If you provide the name of one of the binary resource variables above, this function will
    // return the corresponding original, non-mangled filename (or a null pointer if the name isn't found).
    const char* getNamedResourceOriginalFilename (const char* resourceNameUTF8);
}
