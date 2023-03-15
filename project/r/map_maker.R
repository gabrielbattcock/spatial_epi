map_maker <- function(x) {
  data <- as.data.frame(x)
  data <- data %>%
    mutate(cluster = ifelse(data[,1] > 0 & 
                              data[,2] > 0 & 
                              data[,3] <= 0.05, 
                            "high-high",
                            ifelse(data[,1] <= 0 & 
                                     data[,2] <= 0 & 
                                     data[,3]  <= 0.05, 
                                   "low-low", 
                                   ifelse(data[,1]  > 0 & 
                                            data[,2] <= 0 & 
                                            data[,3] <= 0.05, 
                                          "high-low",
                                          ifelse(data[,1] <= 0 & 
                                                   data[,2] > 0 & 
                                                   data[,3]  <= 0.05,
                                                 "low-high", 
                                                 "p > 0.05")))))
  data <- st_as_sf(data) %>%
    st_transform(4326)
}
