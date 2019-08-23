rm -f $2
echo "//" >> $2
echo "// Do *NOT* edit this file." >> $2
echo "//" >> $2
if which date > /dev/null 2>&1; then
  echo "#define	NGT_BUILD_DATE		\"`date +'%Y/%m/%d %H:%M:%S'`\"" >> $2
fi
if which git > /dev/null 2>&1; then
  echo "#define	NGT_GIT_HASH		\"`git log -1 --format='%H'`\"" >> $2
  echo "#define	NGT_GIT_DATE		\"`git log -1 --format='%cd'`\"" >> $2
  echo "#define	NGT_GIT_TAG		\"`git describe --abbrev=0`\"" >> $2
fi
if which cat > /dev/null 2>&1; then
  echo "#define	NGT_VERSION		\"`cat $1/VERSION`\"" >> $2
fi

touch -r $1/VERSION $2
