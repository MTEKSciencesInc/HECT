Note that some screenshots in the folder are errors/strange behaviour whose exact trigger have yet to be pinpointed.

# Errors/changes
1. Secondary analysis (see screenshot)
- "Error: replacement has length zero" shows up near the bottom of the page
- Having trouble identifying exactly when this can happen

# Possible UI changes 
1. Single Trial Simulation -> Data
- Treatment i labels are not fully visible when there are >= 5

# Quirks
1. Progress bar
- Upon first loading the app, the first run of the Trial design properties' progress bar does not fill properly; it appears and then does not increment until the entire run finishes
- Runs following the first appear to function as expected

2. Type I error (see screenshot)
- "Error: length of 'dimnames' [2] not equal to array extent" sometimes shows up when after running Trial properties using "Compare all arms simultaneously", then changing to "Compare arms against reference treatment" => previous results under Type I error disappear 
- Error disappears once "Run" is pressed again
- This does not always happen consistently