#/usr/bin/env bash
_run_eg_eq_sys_completions()
{
  COMPREPLY=($(compgen -W "$(find $REPO_ROOT/examples -maxdepth 1 -mindepth 1 -type d -exec basename {} \;)" "${COMP_WORDS[1]}"))
}

REPO_ROOT=$( cd -- "$(realpath $( dirname -- "${BASH_SOURCE[0]}" )/..)" &> /dev/null && pwd )
complete -F _run_eg_eq_sys_completions ./scripts/run_eg.sh

EQ_SYS=($(find $REPO_ROOT/examples -maxdepth 1 -mindepth 1 -type d -exec basename {} \;))
for eq_sys in "${EQ_SYS[*]}"
do
    _run_eg_completions()
    {
      COMPREPLY=($(compgen -W "$(find $REPO_ROOT/examples/$eq_sys -maxdepth 1 -mindepth 1 -type d -exec basename {} \;)" "${COMP_WORDS[1]}"))
    }
    complete -F _run_eg_completions "./scripts/run_eg.sh $eq_sys"
done