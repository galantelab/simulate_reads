#
# Rudimentary Bash completion definition for sandy
#

_sandy_database_option() {
	sandy ${@} | awk 'NR>3 && /^\|/ {print $2}'
}

_sandy_database() {
	local cmd="${COMP_WORDS[1]}"
	local subcmd="${COMP_WORDS[2]}"
	local cur="${COMP_WORDS[COMP_CWORD]}"

	local subcmd_opts="
		add
		remove
		restore
		dump
	"

	local long_opts="--help --man"
	local short_opts="-h -u"

	if [[ "$COMP_CWORD" == 2 ]]; then
		case "$cur" in
			--*)
				COMPREPLY=($(compgen -W "$long_opts" -- "$cur"))
				;;
			-*)
				COMPREPLY=($(compgen -W "$short_opts" -- "$cur"))
				;;
			*)
				COMPREPLY=($(compgen -W "$subcmd_opts" -- "$cur"))
				;;
		esac
		return 0
	fi

	long_opts+=" --verbose"
	short_opts+=" -v"

	if [[ "$subcmd" == "add" ]]; then
		case "$cmd" in
			quality)
				long_opts+=" --quality-profile --single-molecule --sequencing-error --source"
				short_opts+=" -s -q -1 -e"
				;;
			expression)
				long_opts+=" --expression-matrix --source"
				short_opts+=" -f -s"
				;;
			variation)
				long_opts+=" --genomic-variation --source"
				short_opts+=" -a -s"
				;;
		esac
	fi

	case "$cur" in
		--*)
			COMPREPLY=($(compgen -W "$long_opts" -- "$cur"))
			;;
		-*)
			COMPREPLY=($(compgen -W "$short_opts" -- "$cur"))
			;;
		*)
			case ${subcmd} in
				add)
					compopt -o default
					COMPREPLY=()
					;;
				remove|dump)
					COMPREPLY=($(compgen -W "$(_sandy_database_option "$cmd")" -- "$cur"))
					;;
				*)
					COMPREPLY=()
					;;
			esac
			;;
	esac
}

_sandy_genome() {
	local cur="${COMP_WORDS[COMP_CWORD]}"
	local prev="${COMP_WORDS[COMP_CWORD-1]}"

	local long_opts="
		--help
		--man
		--verbose
		--prefix
		--output-dir
		--output-format
		--join-paired-ends
		--compression-level
		--append-id
		--id
		--jobs
		--seed
		--coverage
		--sequencing-type
		--quality-profile
		--sequencing-error
		--read-mean
		--read-stdd
		--fragment-mean
		--fragment-stdd
		--genomic-variation
		--genomic-variation-regex
	"

	local short_opts="-h -u -v -p -o -O -1 -x -i -I -j -s -t -q -e -m -d -M -D -A -a -c"

	case "$prev" in
		--genomic-variation|-a)
			COMPREPLY=($(compgen -W "$(_sandy_database_option "variation")" -- "$cur"))
			;;
		--quality-profile|-q)
			COMPREPLY=($(compgen -W "$(_sandy_database_option "quality")" -- "$cur"))
			;;
		--sequencing-type|-t)
			COMPREPLY=($(compgen -W "single-end paired-end" -- "$cur"))
			;;
		--compression-level|-x)
			COMPREPLY=($(compgen -W "$(seq 9)" -- "$cur"))
			;;
		--output-format|-O)
			COMPREPLY=($(compgen -W "bam sam fastq.gz fastq" -- "$cur"))
			;;
		*)
			case "$cur" in
				--*)
					COMPREPLY=($(compgen -W "$long_opts" -- "$cur"))
					;;
				-*)
					COMPREPLY=($(compgen -W "$short_opts" -- "$cur"))
					;;
				*)
					compopt -o default
					COMPREPLY=()
					;;
			esac
			;;
	esac
}

_sandy_transcriptome() {
	local cur="${COMP_WORDS[COMP_CWORD]}"
	local prev="${COMP_WORDS[COMP_CWORD-1]}"

	local long_opts="
		--expression-matrix
		--help
		--man
		--verbose
		--prefix
		--output-dir
		--output-format
		--join-paired-ends
		--compression-level
		--append-id
		--id
		--jobs
		--seed
		--number-of-reads
		--sequencing-type
		--quality-profile
		--sequencing-error
		--read-mean
		--read-stdd
		--fragment-mean
		--fragment-stdd
	"

	local short_opts="-f -h -u -v -p -o -O -1 -x -i -I -j -s -n -t -q -e -m -d -M -D"

	case "$prev" in
		--expression-matrix|-f)
			COMPREPLY=($(compgen -W "$(_sandy_database_option "expression")" -- "$cur"))
			;;
		--quality-profile|-q)
			COMPREPLY=($(compgen -W "$(_sandy_database_option "quality")" -- "$cur"))
			;;
		--sequencing-type|-t)
			COMPREPLY=($(compgen -W "single-end paired-end" -- "$cur"))
			;;
		--compression-level|-x)
			COMPREPLY=($(compgen -W "$(seq 9)" -- "$cur"))
			;;
		--output-format|-O)
			COMPREPLY=($(compgen -W "bam sam fastq.gz fastq" -- "$cur"))
			;;
		*)
			case "$cur" in
				--*)
					COMPREPLY=($(compgen -W "$long_opts" -- "$cur"))
					;;
				-*)
					COMPREPLY=($(compgen -W "$short_opts" -- "$cur"))
					;;
				*)
					compopt -o default
					COMPREPLY=()
					;;
			esac
			;;
	esac
}

_sandy_help() {
	local cmd="${COMP_WORDS[2]}"
	local cur="${COMP_WORDS[COMP_CWORD]}"

	local cmd_opts="
		quality
		expression
		variation
		genome
		transcriptome
	"

	if [[ "$COMP_CWORD" == 2 ]]; then
		COMPREPLY=($(compgen -W "$cmd_opts" -- "$cur"))
		return 0
	fi

	local subcmd_opts="
		add
		remove
		restore
		dump
	"

	case "$cmd" in
		quality|expression|variation)
			COMPREPLY=($(compgen -W "$subcmd_opts" -- "$cur"))
			;;
		*)
			COMPREPLY=()
	esac
}

_sandy() {
	local cmd="${COMP_WORDS[1]}"
	local cur="${COMP_WORDS[COMP_CWORD]}"

	local cmd_opts="
		help
		man
		version
		citation
		quality
		expression
		variation
		genome
		transcriptome
	"

	local long_opts="--help --man"
	local short_opts="-h -u"

	if [[ "$COMP_CWORD" == 1 ]]; then
		case "$cur" in
			--*)
				COMPREPLY=($(compgen -W "$long_opts" -- "$cur"))
				;;
			-*)
				COMPREPLY=($(compgen -W "$short_opts" -- "$cur"))
				;;
			*)
				COMPREPLY=($(compgen -W "$cmd_opts" -- "$cur"))
				;;
		esac
		return 0
	fi

	case "$cmd" in
		help|man)
			_sandy_help
			;;
		quality|expression|variation)
			_sandy_database
			;;
		genome)
			_sandy_genome
			;;
		transcriptome)
			_sandy_transcriptome
			;;
		*)
			COMPREPLY=()
			;;
	esac
}

complete -F _sandy sandy
