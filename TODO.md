# Errors/changes

1. Secondary analysis
- "Error: replacement has length zero" shows up near the bottom of the page

2. Comparison with RCT
- Following a previous run, if settings are changed and Power/Type I error re-estimated, an error is shown under Comparison with RCT ("Error: arguments imply differing number of rows: 5, 2")

3. Single trial design 
- When running with or without "Platform design" checked, sometimes some treatment(s) will only have patients assigned partway through the trial (i.e. not from 0)
- Issue occurs with continuous and binary response
- It's easier to make this behaviour appear with more treatment arms (can set all effect sizes to be the same/trivial)

4. After changing comparison option, results under Power disappear and an error is thrown

# Possible UI changes 
1. Single Trial Simulation -> Data
- Treatment i labels are not fully visible when there are >= 5