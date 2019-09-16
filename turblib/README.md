# Turblib -- A JHU Turbulence Database Cluster C and Fortran Interface Library

See the end of the file for license conditions.

## About this package

This library provides a C and Fortran Interface wrapper around gSOAP for calling the [JHU Turbulence Database Cluster](http://turbulence.pha.jhu.edu/).
More information can be found at: [http://turbulence.pha.jhu.edu/help/c-fortran/](http://turbulence.pha.jhu.edu/help/c-fortran/)

`make` will build both the Fortran and C sample code. Use `make DEMO_turbc` if you do not have a Fortran compiler installed.  Note: On some platforms, you may need to use *gmake* (GNU Make) instead.

## Identification Token (required for large queries)
While our service is open to anyone, we would like to keep track of who is using the service, and how. To this end, we would like each user or site to obtain an authorization token from us: [http://turbulence.pha.jhu.edu/help/authtoken.html](http://turbulence.pha.jhu.edu/help/authtoken.html)

If you are just experimenting, the default token included in these test files should be valid.

## Error Handling

By default, the turbulence library will print errors out to the console and abort the program. You can change this behavior by calling `turblibSetExitOnError(0)` before making any calls to the turbulence database. You will then be responsible for checking the return codes.

#### C Example:
```
turblibSetExitOnError(0);
...
if (getVelocity ( ... ) != SOAP_OK) {
  turblibPrintError();
  exit(1);
}
```

#### Fortran Example:
```
integer, parameter :: SOAP_OK = 0  ! From stdsoap2.h
integer rc
CALL turblibSetExitOnError(0)
...
rc = getVelocity( ... )
if (rc.ne.SOAP_OK) then
  CALL turblibPrintError()
  STOP
end if
```

## Handling intermittent failures

We suggest writing long-running code to handle any intermittent network or database failures which may occur.  The following sample code will try each query up to 30 times, with a delay of 1 minute between each attempt.

In addition, we suggest adding checkpoints to your code to allow for the continuation after longer failures.

#### C:
```
int attempts = 0;
turblibSetExitOnError(0);
while (getVelocity ( ... ) != SOAP_OK) {
  if (attempts++ > 30) {
    printf("Fatal Error: too many failures\n");
    exit(1);
  } else {
    printf("Temporary Error: %s\n", turblibGetErrorString());
  }
  sleep(60);
}
```

#### Fortran:
```
integer attempts
attempts = 0
CALL turblibSetExitOnError(0)
do while (getVelocity( ... ).ne.0)
  attempts = attempts + 1
  if (attempts.ge.30) then
    write(*,*) 'Fatal error: too many failures'
    CALL turblibPrintError()
    STOP
  else
    write(*,*) 'Temporary Error (#', attempts, '):'
    CALL turblibPrintError()
    CALL sleep(1)
  end if
end do
```

## License
Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at [http://www.apache.org/licenses/LICENSE-2.0](http://www.apache.org/licenses/LICENSE-2.0)

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.