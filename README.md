# RelaxiFI
RelaxiFI is an open-source python3 package designed to model the stretching of CO2-dominated fluid inclusions by stress-relaxation dislocation creep, it is based on the model and parameters described by Wanamaker and Evans (1989). Currently, RelaxiFI is not equipped to model stretching in minerals other than olivine, sorry!

NOTE: Coolprop (http://www.coolprop.org/) available at PyPI (https://pypi.org/project/CoolProp/) is required to perform calculations using the EOS of Span and Wagner (1996) in this package. Make sure to install it using pip install CoolProp or !pip install CoolProp (if from a jupyter notebook)

For more information and usage examples, please see https://relaxifi.readthedocs.io/

Copyright (C) <2023 Charlotte Devitre> 

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version. 
You cannot use RelaxiFI to create closed source software.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details (https://github.com/cljdevitre/RelaxiFI/blob/main/LICENSE) or see <https://www.gnu.org/licenses/>.