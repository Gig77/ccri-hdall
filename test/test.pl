use strict;
use warnings FATAL => qw( all );

use lib "/mnt/projects/generic/scripts";
use Generic;
use Log::Log4perl qw(:easy);

DEBUG("debug");
INFO("info");
WARN("warn");
ERROR("error");
FATAL("fatal");
