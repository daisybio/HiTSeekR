positionsToAlphaName <- function(positions){
  if(max(positions) == 96)
  {
    alpha <- alphaNames(row=8, column=12)
  }
  else if(max(positions) == 384)
  {
    alpha <- alphaNames(row=16, column=24)
  }
  else if(max(positions) == 1536)
  {
    alpha <- alphaNames(row=32, column=48)
  }
  
  return(alpha[positions])
}