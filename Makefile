run:
	g++ -DDEBUG="true" -DTYPES="FAST_FIXED(20, 10), FIXED(20, 10),FIXED(20, 20),FLOAT" -DSIZES="S(1920,1080),S(36,84)" -std=c++2b fluid.cpp
	./a.out --p-type="FLOAT" --v-type="FAST_FIXED(20, 10)" --v-flow-type="FIXED(20, 10)" --n=36 --m=84 --field-file="field"

runSave:
	-mkdir save
	g++ -DDEBUG="true" -DTYPES="FAST_FIXED(20, 10), FIXED(20, 10),FIXED(20, 20),FLOAT" -DSIZES="S(1920,1080),S(36,84)" -DSAVE -std=c++2b fluid.cpp
	./a.out --p-type="FLOAT" --v-type="FAST_FIXED(20, 10)" --v-flow-type="FIXED(20, 10)" --n=36 --m=84 --field-file="field"

runLoad:
	g++ -DDEBUG="true" -DTYPES="FAST_FIXED(20, 10), FIXED(20, 10),FIXED(20, 20),FLOAT" -DSIZES="S(1920,1080),S(36,84)" -std=c++2b fluid.cpp
	./a.out --p-type="FLOAT" --v-type="FAST_FIXED(20, 10)" --v-flow-type="FIXED(20, 10)" --n=36 --m=84 --field-file="field" --load=165

testBaseCode:
	g++ -DDEBUG="true" -std=c++2b base_fluid.cpp
	./a.out

CreateBuildLib:
	-mkdir fluid_build
	cd fluid_build && cmake ..

BuildRunCmakeRelease:
	cd fluid_build && cmake .. -DSetDebugBuildType=OFF -DTYPES="FAST_FIXED(20, 10), FIXED(20, 10),FIXED(20, 20),FLOAT" -DSIZES="S(1920,1080),S(36,84)"
	cd fluid_build && cmake --build .
	cd fluid_build && ./Fluid3Project --p-type="FAST_FIXED(20, 10)" --v-type="FAST_FIXED(20, 10)" --v-flow-type="FAST_FIXED(20, 10)" --n=36 --m=84 --field-file="../field"

BuildRunCmakeReleaseSemiLarge:
	cd fluid_build && cmake .. -DSetDebugBuildType=OFF -DTYPES="FAST_FIXED(20, 10)" -DSIZES="S(60, 150)"
	cd fluid_build && cmake --build .
	cd fluid_build && ./Fluid3Project --p-type="FAST_FIXED(20, 10)" --v-type="FAST_FIXED(20, 10)" --v-flow-type="FAST_FIXED(20, 10)" --n=60 --m=150 --field-file="../field_semilarge"

BuildRunCmakeReleaseLarge:
	cd fluid_build && cmake .. -DSetDebugBuildType=OFF -DTYPES="FAST_FIXED(20, 10)"
	cd fluid_build && cmake --build .
	cd fluid_build && ./Fluid3Project --p-type="FAST_FIXED(20, 10)" --v-type="FAST_FIXED(20, 10)" --v-flow-type="FAST_FIXED(20, 10)" --n=64 --m=512 --end-on-tick=200 --field-file="../field_large"

BuildRunCmakeDebug:
	cd fluid_build && cmake -S .. -DSetDebugBuildType=ON -DTYPES="FAST_FIXED(20, 10), FIXED(20, 10),FIXED(20, 20),FLOAT" -DSIZES="S(1920,1080),S(36,84)"
	cd fluid_build && cmake --build .

	cd fluid_build && ./Fluid3Project --p-type="FAST_FIXED(20, 10)" --v-type="FAST_FIXED(20, 10)" --v-flow-type="FAST_FIXED(20, 10)" --n=36 --m=84 --field-file="../field"
