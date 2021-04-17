//macro to evaluate density of Ki67 positive cells (cells/mm2)

//open a file that has to be analysed and choose the series

run("Set Measurements...", "area mean integrated stack limit display nan redirect=None decimal=2");

titles = getList("image.titles");

for(i=0; i<lengthOf(titles); i++) {

selectWindow(titles[i]);
print(titles[i]);
//processing on DAPI channel
Stack.setChannel(2);
run("Duplicate...", "duplicate channels=2");
run("Subtract Background...", "rolling=200 slice");
run("Gaussian Blur...", "sigma=8 slice");

setAutoThreshold("Huang dark");
waitForUser("Check threshold on DAPI");

getThreshold(lower, upper);
print("DAPI low thresh:"+lower);

//run("Convert to Mask", "method=Default background=Dark only black");
setOption("BlackBackground", true);
run("Convert to Mask");

//waitForUser("");
run("Measure");

rename("Mask DAPI");

//process on KI67 channel

selectWindow(titles[i]);
Stack.setChannel(1);

run("Duplicate...", "duplicate channels=1");
run("Subtract Background...", "rolling=50 slice");
run("Gaussian Blur...", "sigma=2 slice");

setAutoThreshold("Otsu dark");
waitForUser("Check threshold on KI67");

getThreshold(lower, upper);
print("KI67 low thresh:"+lower);

//run("Convert to Mask", "method=Triangle background=Dark only black");
setOption("BlackBackground", true);
run("Convert to Mask");

run("Open");
run("Close-");
//run("Watershed");

kimask=getTitle();

imageCalculator("AND", "Mask DAPI", kimask);
run("Convert to Mask");
run("Measure");
//waitForUser("");

//run("Analyze Particles...", "size=25-200 exclude include summarize");
selectWindow("Mask DAPI");
run("Close");
selectWindow(kimask);
run("Close");

DAPI_area = getResult("Area", 2*i);
KI67_area = getResult("Area", 2*i+1);

pos_f = KI67_area/DAPI_area;
print("Area KI67:" + KI67_area);
print("Area DAPI:" + DAPI_area);
print("Area ratio KI67/DAPI:" + pos_f);

//waitForUser("Check results Pause");
}
