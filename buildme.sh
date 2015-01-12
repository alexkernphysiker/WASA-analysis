#/bin/bash
echo "COMPILING..."
(cd Analysis_test;qmake Analysis_test.pro;make clean; make)
if [ $# -ne 1 ]; then
    echo "Binaries are not copied"
    exit
fi
echo "Copying the binaries to " $1
find -name *.so -exec cp '{}' $1 \;
find -name *.so.* -exec cp '{}' $1 \;
find -name *_app -exec cp '{}' $1 \;
find -name *_app.* -exec cp '{}' $1 \;
echo "done."
echo "Before running programs you may need to execute:"
echo "export LD_LIBRARY_PATH=."
