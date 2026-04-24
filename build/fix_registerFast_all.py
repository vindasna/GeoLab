import re

# 1. Update bundlesData.h
with open("../src/bundlesData.h", "r") as f:
    h_content = f.read()

enum_def = """enum class TransformType {
    RIGID,
    AFFINE
};

struct Cluster {"""

# Insert enum before struct Cluster if it's not there
if "enum class TransformType" not in h_content:
    h_content = h_content.replace("struct Cluster {", enum_def)

# Replace old registerFast declaration
old_declaration = "void registerFast(BundlesData& target, float qbThreshold = 15.0f, int iterations = 10, float outlierRejectionPct = 0.20f);"
new_declaration = """Eigen::Matrix4f registerFast(BundlesData& target, 
                                  TransformType regType = TransformType::RIGID, 
                                  float qbThreshold = 10.0f, 
                                  int iterations = 10, 
                                  float outlierRejectionPct = 0.20f,
                                  BundlesData* outMovingCentroids = nullptr,
                                  BundlesData* outTargetCentroids = nullptr);"""

# The existing declaration in .h might be different based on previous command output. 
# Let's see what was exactly in bundlesData.h for registerFast.
# It was: void registerFast(BundlesData& target, float qbThreshold = 15.0f, int iterations = 10, float outlierRejectionPct = 0.20f);
h_content = re.sub(
    r'void\s+registerFast\s*\([^)]+\)\s*;',
    new_declaration,
    h_content
)

with open("../src/bundlesData.h", "w") as f:
    f.write(h_content)


# 2. Update bundlesData.cxx
with open("../src/bundlesData.cxx", "r") as f:
    cxx_content = f.read()

# Remove enum class TransformType from .cxx
cxx_enum = """// Enum for the user to select the registration type
enum class TransformType {
    RIGID,
    AFFINE
};"""

if cxx_enum in cxx_content:
    cxx_content = cxx_content.replace(cxx_enum, "")
else:
    cxx_content = re.sub(
        r'enum class TransformType\s*\{[^}]+\}\s*;',
        "",
        cxx_content
    )

with open("../src/bundlesData.cxx", "w") as f:
    f.write(cxx_content)


# 3. Update apps/registerTractogramsFast.cxx
with open("../apps/registerTractogramsFast.cxx", "r") as f:
    app_content = f.read()

# Add -regType argument parsing
argument_position = """    index_output = getFlagPosition( argc, argv, "-o" ) ;
    index_qbThr = getFlagPosition( argc, argv, "-qbThr" ) ;"""

new_argument = """    index_output = getFlagPosition( argc, argv, "-o" ) ;
    int index_regType = getFlagPosition( argc, argv, "-regType" ) ;
    index_qbThr = getFlagPosition( argc, argv, "-qbThr" ) ;"""

if argument_position in app_content:
    app_content = app_content.replace(argument_position, new_argument)

# Add help message
help_message = """                  << "-o : Output registered moved tractogram path\\n"
                  << "[-qbThr] : Distance threshold (mm) for QuickBundles centroids (default = 10.0)\\n\""""

new_help = """                  << "-o : Output registered moved tractogram path\\n"
                  << "[-regType] : Type of registration (rigid or affine, default = rigid)\\n"
                  << "[-qbThr] : Distance threshold (mm) for QuickBundles centroids (default = 10.0)\\n\""""

if help_message in app_content:
    app_content = app_content.replace(help_message, new_help)

# Add variable assignment
var_assign = """    float qbThr = 10.0f ;"""
new_var_assign = """    TransformType regType = TransformType::RIGID ;
    float qbThr = 10.0f ;"""

if var_assign in app_content:
    app_content = app_content.replace(var_assign, new_var_assign)

parse_logic = """    if ( index_qbThr ) qbThr = std::stof( argv[ index_qbThr + 1 ] ) ;"""
new_parse_logic = """    if ( index_regType ) {
        std::string typeStr = argv[ index_regType + 1 ] ;
        if ( typeStr == "affine" || typeStr == "Affine" ) regType = TransformType::AFFINE ;
        else if ( typeStr == "rigid" || typeStr == "Rigid" ) regType = TransformType::RIGID ;
        else { std::cout << "Invalid -regType argument: " << typeStr << ", using default Rigid" << std::endl; }
    }
    if ( index_qbThr ) qbThr = std::stof( argv[ index_qbThr + 1 ] ) ;"""

if parse_logic in app_content:
    app_content = app_content.replace(parse_logic, new_parse_logic)

# Print verbose
print_verbose = """        std::cout << "Output Path      : " << outputPath << std::endl ;"""
new_print_verbose = """        std::cout << "Output Path      : " << outputPath << std::endl ;
        std::cout << "Registration Type: " << (regType == TransformType::RIGID ? "Rigid" : "Affine") << std::endl ;"""

if print_verbose in app_content:
    app_content = app_content.replace(print_verbose, new_print_verbose)

# Function call update
call_update = """    movingBundles.registerFast( fixedBundles, qbThr, nbIter, outlierRejectionThr ) ;"""
new_call_update = """    Eigen::Matrix4f transformMatrix = movingBundles.registerFast( fixedBundles, regType, qbThr, nbIter, outlierRejectionThr ) ;
    if (verbose) {
        std::cout << "Cumulative Transformation Matrix:\\n" << transformMatrix << std::endl ;
    }"""

if call_update in app_content:
    app_content = app_content.replace(call_update, new_call_update)

with open("../apps/registerTractogramsFast.cxx", "w") as f:
    f.write(app_content)

