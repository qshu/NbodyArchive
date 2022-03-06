# The purpose of this script is to
# generate the job description.
# write_rsl() is called in submit.sh

function write_rsl() {
cat > $1 <<EOF
<?xml version="1.0" encoding="UTF-8"?>

<job>
  <executable>nb6_wrapper_${task_id}.sh</executable>
  <!-- <directory>\${GLOBUS_USER_HOME}</directory> -->
  
  <stdout>nb6_wrapper_stdout_${task_id}</stdout>
  <stderr>nb6_wrapper_stderr_${task_id}</stderr>
EOF

if [ x$JOB_QUEUE != 'xnot_specified' ]; then

cat >> $1 <<EOF
  <queue>$JOB_QUEUE</queue>
EOF
fi

cat >> $1 <<EOF
  <fileStageIn>

  <transfer>
  <sourceUrl>${THIS_URL}/tmp/nb6_wrapper.sh_${task_id}</sourceUrl>
  <destinationUrl>file:///\${GLOBUS_USER_HOME}/nb6_wrapper_${task_id}.sh</destinationUrl>
  </transfer>

  <transfer>
  <sourceUrl>${THIS_URL}/tmp/nbody6_src_${task_id}.tar.gz</sourceUrl>
  <destinationUrl>file:///\${GLOBUS_USER_HOME}/nbody6_src_${task_id}.tar.gz</destinationUrl>
  </transfer>

  <transfer>
    <sourceUrl>${THIS_URL}/libexec/hosts.env</sourceUrl>
    <destinationUrl>file:///\${GLOBUS_USER_HOME}/hosts.env_${task_id}</destinationUrl>
  </transfer>

EOF

if [ x$RESTART = 'xyes' ]; then

cat >> $1 <<EOF
  <transfer>
    <sourceUrl>${THIS_URL}/${COMMON_BLOCK_FILE}</sourceUrl>
    <destinationUrl>file:///\${GLOBUS_USER_HOME}/comm.1_${task_id}</destinationUrl>
  </transfer>

EOF
fi

if [ x$HAVE_INITIAL_DATA_FILE = 'xyes' ]; then

cat >> $1 <<EOF
  <transfer>
    <sourceUrl>${THIS_URL}/${INITIAL_DATA_FILE}</sourceUrl>
    <destinationUrl>file:///\${GLOBUS_USER_HOME}/dat.10_${task_id}</destinationUrl>
  </transfer>

EOF
fi


if [ -f libexec/plugins/* ]; then
plugin_no=1
for i in libexec/plugins/*; do
cat >> $1 <<EOF
  <transfer>
    <sourceUrl>${THIS_URL}/${i}</sourceUrl>
    <destinationUrl>file:///\${GLOBUS_USER_HOME}/nb6_plugin_${plugin_no}_${task_id}</destinationUrl>
  </transfer>

EOF
plugin_no=$((plugin_no+1))
done
fi

cat >> $1 <<EOF
  <transfer>
  <sourceUrl>${THIS_URL}/${PARAMETER_INPUT_FILE}</sourceUrl>
  <destinationUrl>file:///\${GLOBUS_USER_HOME}/parameter.in_${task_id}</destinationUrl>
  </transfer>
  </fileStageIn>

  <fileStageOut>
  <transfer>
  <sourceUrl>file:///\${GLOBUS_USER_HOME}/nb6_wrapper_stdout_${task_id}</sourceUrl>
  <destinationUrl>${THIS_URL}/outfiles/${task_id}/wrapper.out</destinationUrl>
  </transfer>
 
  <transfer>
    <sourceUrl>file:///\${GLOBUS_USER_HOME}/nb6_wrapper_stderr_${task_id}</sourceUrl>
    <destinationUrl>${THIS_URL}/outfiles/${task_id}/wrapper.err</destinationUrl>
  </transfer>

  <transfer>
    <sourceUrl>file:///\${GLOBUS_USER_HOME}/nbody6.out_${task_id}</sourceUrl>
    <destinationUrl>${THIS_URL}/outfiles/${task_id}/nbody6.out</destinationUrl>
  </transfer>
  
  <transfer>
    <sourceUrl>file:///\${GLOBUS_USER_HOME}/comm.1_${task_id}</sourceUrl>
    <destinationUrl>${THIS_URL}/outfiles/${task_id}/comm.1</destinationUrl>
  </transfer>

   <transfer>
    <sourceUrl>file:///\${GLOBUS_USER_HOME}/comm.2_${task_id}</sourceUrl>
    <destinationUrl>${THIS_URL}/outfiles/${task_id}/comm.2</destinationUrl>
  </transfer>

   <transfer>
    <sourceUrl>file:///\${GLOBUS_USER_HOME}/conf.3_${task_id}</sourceUrl>
    <destinationUrl>${THIS_URL}/outfiles/${task_id}/conf.3</destinationUrl>
  </transfer>
  
   <transfer>
    <sourceUrl>file:///\${GLOBUS_USER_HOME}/hia.12_${task_id}</sourceUrl>
    <destinationUrl>${THIS_URL}/outfiles/${task_id}/hia.12</destinationUrl>
  </transfer>

   <transfer>
    <sourceUrl>file:///\${GLOBUS_USER_HOME}/lagr.7_${task_id}</sourceUrl>
    <destinationUrl>${THIS_URL}/outfiles/${task_id}/lagr.7</destinationUrl>
  </transfer>

  </fileStageOut>

  <fileCleanUp>
  
  <deletion>
    <file>file:///\${GLOBUS_USER_HOME}/nb6_wrapper_${task_id}.sh</file>
  </deletion>
 
  <deletion>
    <file>file:///\${GLOBUS_USER_HOME}/nbody6_src_${task_id}.tar.gz</file>
  </deletion>

  <deletion>
    <file>file:///\${GLOBUS_USER_HOME}/hosts.env_${task_id}</file>
  </deletion>
  
  <deletion>
    <file>file:///\${GLOBUS_USER_HOME}/parameter.in_${task_id}</file>
  </deletion>

  <deletion>
    <file>file:///\${GLOBUS_USER_HOME}/nb6_wrapper_stdout_${task_id}</file>
  </deletion>

  <deletion>
    <file>file:///\${GLOBUS_USER_HOME}/nb6_wrapper_stderr_${task_id}</file>
  </deletion>

  <deletion>
    <file>file:///\${GLOBUS_USER_HOME}/nbody6.out_${task_id}</file>
  </deletion>

  <deletion>
    <file>file:///\${GLOBUS_USER_HOME}/comm.1_${task_id}</file>
  </deletion>

  <deletion>
    <file>file:///\${GLOBUS_USER_HOME}/comm.2_${task_id}</file>
  </deletion>

  <deletion>
    <file>file:///\${GLOBUS_USER_HOME}/conf.3_${task_id}</file>
  </deletion>

  <deletion>
    <file>file:///\${GLOBUS_USER_HOME}/hia.12_${task_id}</file>
  </deletion>

  <deletion>
    <file>file:///\${GLOBUS_USER_HOME}/lagr.7_${task_id}</file>
  </deletion>

EOF
  
  plugin_no=1
  for i in libexec/plugins/*; do
    if [ ${i} == 'libexec/plugins/*' ]; then
        break
     else    
  cat >> $1 <<EOF
    <deletion>
        <file>file:///\${GLOBUS_USER_HOME}/nb6_plugin_${plugin_no}_${task_id}</file>
    </deletion>

EOF
  plugin_no=$((plugin_no+1))
    fi
  done

cat >> $1 <<EOF
  </fileCleanUp>

<!--
<extensions>
<gw>
HOSTNAME="gavo2.aip.de"
</gw>
</extensions>
-->

</job>
EOF
}

