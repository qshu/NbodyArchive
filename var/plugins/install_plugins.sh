#!/bin/sh

if [ -f "../../scripts/common.sh" ]; then
  . ../../scripts/common.sh
else
  echo "Error. Could not find common.sh"
  exit 1
fi

PLUGINS=demo_plugin


function copy_plugin()
{
  if [ -f $1 ];then
    CP="cp $1 ../../../libexec/plugins"
    echo "$CP"
    $CP
    if [ $? -eq 0 ]; then
      echo "Plugin $(basename $1) sucessfully installed."
    else
      echo "Could not copy $(basename $1) to plugin directory."
      exit 1
    fi
  else
    echo "install_plugin(): can not find $1"
    exit 1 
  fi
}



function install_plugins()
{
# Replace the #SOURCE_LINE_MARKER in the plugins
# with the real line number
cd src
for i in $@; do
  if [ ! -f $i ]; then
    echo "No such plugin: $i"
    continue
  fi

  echo ""
  echo "--------------------------"
  echo "Processing $i"
  echo "--------------------------"

  search_for_marker=`grep -n '#SOURCE_LINE_MARKER' $i > line.tmp`
  if [ $? -ne 0 ]; then
    echo "No #SOURCE_LINE_MARKER found in $i"
    copy_plugin $i
    continue
  fi
  get_line=`awk -F: '{print $1}' line.tmp`
  
  line=$((get_line + 1))
  echo "Found #SOURCE_LINE_MARKER in $i at $get_line"
  sed -e "s/\#SET_SOURCE_LINE/$line/" $i > ../$i.bin
  [ -f ../$i.bin ] && chmod +x ../$i.bin
  echo "Plugin ${i}.bin generated."
  copy_plugin "../${i}.bin"
  [ -f ../${i}.bin ] && rm ../${i}.bin
done

[ -f line.tmp ] && rm line.tmp
}


#if [ $# -lt 1 ]; then
#  echo "usage: $0 plugin1 plugin2...plugin[n]"
#fi

echo "Which plugin(s) should be installed/updated?"

echo "--------------------"

for i in src/*; do
echo $(basename ${i})
done

echo "--------------------"
echo ""

read -p "Enter: " plugs

install_plugins $plugs


