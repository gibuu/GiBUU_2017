program test_XsectionRatios

use inputGeneral, only : readInputGeneral
use XsectionRatios


implicit none

call readInputGeneral
call init_Database

call init

end program test_XsectionRatios
