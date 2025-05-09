#!/bin/sh
skip=49

tab='	'
nl='
'
IFS=" $tab$nl"

umask=`umask`
umask 77

gztmpdir=
trap 'res=$?
  test -n "$gztmpdir" && rm -fr "$gztmpdir"
  (exit $res); exit $res
' 0 1 2 3 5 10 13 15

case $TMPDIR in
  / | /*/) ;;
  /*) TMPDIR=$TMPDIR/;;
  *) TMPDIR=/tmp/;;
esac
if type mktemp >/dev/null 2>&1; then
  gztmpdir=`mktemp -d "${TMPDIR}gztmpXXXXXXXXX"`
else
  gztmpdir=${TMPDIR}gztmp$$; mkdir $gztmpdir
fi || { (exit 127); exit 127; }

gztmp=$gztmpdir/$0
case $0 in
-* | */*'
') mkdir -p "$gztmp" && rm -r "$gztmp";;
*/*) gztmp=$gztmpdir/`basename "$0"`;;
esac || { (exit 127); exit 127; }

case `printf 'X\n' | tail -n +1 2>/dev/null` in
X) tail_n=-n;;
*) tail_n=;;
esac
if tail $tail_n +$skip <"$0" | gzip -cd > "$gztmp"; then
  umask $umask
  chmod 700 "$gztmp"
  (sleep 5; rm -fr "$gztmpdir") 2>/dev/null &
  "$gztmp" ${1+"$@"}; res=$?
else
  printf >&2 '%s\n' "Cannot decompress $0"
  (exit 127); res=127
fi; exit $res
�s��gtest.sh �T[KA~�_q�lĐn.B����P���H�d'�-��twRi�FmJ���jQD*i�DKo���\֧����f5���,g��o�sf|}ጪ�3�� ���1>5
�0u�Btl��b n-V�K���Fk�v��U
Q�S�Ǆ�Tb�ߍa��`��}�r\[hmU��k�����Ls����kt�4�.W?5�_����ڌ����q���ު��̷g�[�~ڻ����.� O�x�AV����pBoo|�l6:��;��7�Nu��X�^����H���)i���tX�4,���J��Djrl)F����x�A�QXp�c�̃�������tW�Q���z�;� b��{
W(�h���	�� �z(��p~AL2&�H�B���H�H2��c�`0E艑Q�0��*etC&��9�o�,%��C�J��J)��e´J����a2�Z F��nE"D4�h9A�N��L�b7�VT�0���8������283(t����^����)��Uq�<� x�@ ɆN���Rɹ>Nte0/�
�o�e�Q����'t��l�wZ�ˎ�`ϭ�׶;�s��w�Z��\k�o���q���ve^O$�V,Ѵ����%�^(�J�h��'�?��N^T�1#g�tψ�R8a��yi&9�Įmn&�
��Q({�6vgb�Q��gy/�']������j ��ʖL��4�'��e�A�$��_�#����O1O�f�����x�H�y���b��=����ݗ�J�#�Uݻ)^s��`T��!F�	������  