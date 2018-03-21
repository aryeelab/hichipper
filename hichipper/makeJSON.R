library(RJSONIO)

s1 <- data.frame(
  type="longrange",
  url=URL_SAMPLE1,
  name=SAMPLE1,
  mode="arc",
  colorpositive="#00441B",  # change colors as appropriate
  height="50"
)

s2 <- data.frame(
  type="longrange",
  url=URL_SAMPLE2,
  name=SAMPLE2,
  mode="arc",
  colorpositive="#00441B", # change colors as appropriate
  height="50"
)

exportJson <- toJSON(do.call(rbind, list(s1,s2)))
write(exportJson, "washu_tracks.json")
