
setClass("CANCERTable",
         slots=list(
           count="data.frame",
           metadata="data.frame",
           gene_type="vector"
         ),
         prototype=list(
           count=data.frame(matrix(nr=0,nc=0)),
           metadata=data.frame(matrix(nr=0,nc=0)),
          gene_type=vector()
         )
)

setClass("lnRAN_mRNA",
         slots=list(
           lnRNA="data.frame",
           mRNA="data.frame",
           metadata="data.frame"
         ),
         prototype=list(
           lnRNA=data.frame(matrix(nr=0,nc=0)),
           mRNA=data.frame(matrix(nr=0,nc=0)),
           metadata=data.frame(matrix(nr=0,nc=0))
         )
)
setClass("DEGTable",
         slots=list(
           deg="data.frame",
           gene_type="vector"
         ),
         prototype=list(
           deg=data.frame(matrix(nr=0,nc=0)),
           gene_type=vector()
         )
)



setClass("CoxMulti",
         slots=list(
           single="data.frame",
           multi="data.frame",
           data="data.frame"
         ),
         prototype=list(
           single=data.frame(matrix(nr=0,nc=0)),
           multi=data.frame(matrix(nr=0,nc=0)),
           data=data.frame(matrix(nr=0,nc=0))
         )
)



