# bash completion for dkms
# Copied from the Mandriva dkms package

# This function complete on available kernels
#
#_kernels()
#{
#        COMPREPLY=( $( cd /lib/modules && compgen -d -- "$cur" ) )
#}

# complete on full directory names under $1
#_subdirectories()
#{
#	COMPREPLY=( $( cd $1 && compgen -d -- "$cur" ) )
#}

# complete on $2 part of filenames matching pattern $1 under /usr/src
#_filename_parts()
#{
#	COMPREPLY=( $( runMode ls -F /usr/src/ 2>/dev/null | grep -E '^'$1'/$' \
#		| sed -r -e 's/^([^-]+)-(.+)\/$/\'$2'/' | grep "^$cur" ) )
#}

_3seq()
{
	local cur prev runMode module i

	COMPREPLY=()
	cur=${COMP_WORDS[COMP_CWORD]}

	if [[ $COMP_CWORD -eq 1 ]] ; then
		COMPREPLY=( $( compgen -W "-help -version -gen-p -check-p -write -info \
                                   -match -triplet -single -full" -- $cur ) )

#	elif [[ $COMP_CWORD -eq 2 ]] ; then
#        case $runMode in
#				add)
#					options='-c --rpm_safe_upgrade'
#					;;
#				remove)
#					options='--rpm_safe_upgrade'
#					;;
#				build)
#					options='--config'
#					;;
#				mkdriverdisk)
#					options='-d --distro -r --release --size'
#					;;
#				ldtarball)
#					options='--archive --force'
#					;;
#				mktarball)
#					options='--source-only --binaries-only'
#					;;
#				mkrpm)
#					options='--source-only'
#					;;
#				mkkmp)
#					options='--spec'
#					;;
#				match)
#					options='--templatekernel'
#    				;;
#			esac
#
#            options="$options -m -v -k -a --arch -q --quiet -V \
#				--version --all --no-prepare-kernel \
#				--no-clean-kernel --kernelsourcedir \
#				--directive"
#
#			COMPREPLY=( $( compgen -W "$options" -- $cur ) )
#
#    else
#
#		prev=${COMP_WORDS[COMP_CWORD-1]}
#		runMode=${COMP_WORDS[1]}
#		case $prev in
#			-m)
#				if [ "$runMode" = 'add' ]; then
#					_filename_parts '.*-.*' 1
#				else
#					_subdirectories /var/lib/dkms
#				fi
#				return 0
#				;;
#			-v)
#				for (( i=1; i < COMP_CWORD; i++ )); do
#					if [[ "${COMP_WORDS[i]}" == -m ]]; then
#						module=${COMP_WORDS[i+1]}
#						break
#					fi
#				done
#				if [ -n "$module" ]; then
#					if [ "$runMode" = 'add' ]; then
#						_filename_parts "$module-.*" 2
#					else
#						_subdirectories /var/lib/dkms/$module
#					fi
#					return 0
#				fi
#				;;
#			-k)
#				_kernels
#				return 0
#				;;
#			-@\(c|-spec|-archive|-config\))
#				_filedir
#				return 0
#				;;
#			--kernelsourcedir)
#				_filedir -d
#				return 0
#				;;
#		esac
#
#
#		if [[ "$cur" == -* ]]; then
#			case $runMode in
#				add)
#					options='-c --rpm_safe_upgrade'
#					;;
#				remove)
#					options='--rpm_safe_upgrade'
#					;;
#				build)
#					options='--config'
#					;;
#				mkdriverdisk)
#					options='-d --distro -r --release --size'
#					;;
#				ldtarball)
#					options='--archive --force'
#					;;
#				mktarball)
#					options='--source-only --binaries-only'
#					;;
#				mkrpm)
#					options='--source-only'
#					;;
#				mkkmp)
#					options='--spec'
#					;;
#				match)
#					options='--templatekernel'
#					;;
#			esac
#
#			options="$options -m -v -k -a --arch -q --quiet -V \
#				--version --all --no-prepare-kernel \
#				--no-clean-kernel --kernelsourcedir \
#				--directive"
#
#			COMPREPLY=( $( compgen -W "$options" -- $cur ) )
#		fi
	fi
}
complete -F _3seq 3seq