# Este script está diseñado para tomar la base de datos de Distancias.db, separar todos los datos
# para metilaciones que estén dentro del gen según caigan en exones o intrones. Finalmente, a falta de una
# métrica mejor, sacamos el promedio para cada uno. Así, cada gen tiene un score para la metilación
# que se encuentra antes, en sus exones, en sus intrones y después de él, correspondiente al promedio de
# los scores de las regiones metiladas que pertenecen a cada categoría.....

library(RSQLite)
library(argparse)

parser = ArgumentParser(description = 'Get database')
parser$add_argument('database', help='the database')

args = parser$parse_args()

# nos conectamos a la base de datos que nos interesa
drv = dbDriver('SQLite')
db = dbConnect(drv, args$database)

# nos conectamos a la base de datos que tiene la información sobre los exones
drv2 = dbDriver('SQLite')
db2 = dbConnect(drv2, 'Homo_sapiens.GRCh38.80.db')

# nos conectamos a la base de datos que tiene la información de expresión
drv3 = dbDriver('SQLite')
db3 = dbConnect(drv3, 'Modifications_RPKM.db')

#############################################################################################

separar_ex_int = function(tabla){
    # cargamos los datos, ignorando los NA
    datos = dbGetQuery(db, paste('select * from ', tabla, ' where closest_gene_name!="NA" and block_location="inside"', sep =''))
    
    # guardamos el tejido del que estamos hablando
    # tomamos los primeros cuatro caracteres del nombre de la tabla (id del tejido)
    tejido = substr(tabla, 1, 4)
    # reasignamos la variable tejido con el nombre común de acuerdo al id
    if(tejido == 'E066'){
        tejido = 'Liver'
    }else if(tejido == 'E096'){
        tejido = 'Lung'
    }else if(tejido == 'E097'){
        tejido = 'Ovary'
    }else if(tejido == 'E098'){
        tejido = 'Pancreas'
    }
    
    # inicializamos los archivos donde guardaremos los exones e intrones antes de cargarlos a la base de datos
    tmpfile_exones = file(paste(tejido, '_avg_exones.txt', sep=''))
    tmpfile_intrones = file(paste(tejido, '_avg_intrones.txt', sep=''))
    

    # juntamos la lista de genes
    genes = names(table(datos$closest_gene_name))
    
    # inicializamos el número de errores en cero
    num_errores = 0
 
    # iteramos por cada gen
    for(gen in genes){
        # tomamos los datos para ese gen
        dummy = dbGetQuery(db, paste('select * from ', tabla, ' where closest_gene_name="', gen, '" and block_location="inside"', sep = ''))
        
        # guardamos el cromosoma de ese gen
        chr = dummy[1,'chr']
        
        # reinicializamos las variables para exones e intrones
        exones = c()
        intrones = c()
        
        # iteramos por cada bloque para ese gen
        for(line in 1:dim(dummy)[1]){
            # guardamos el inicio y el final del bloque
            inicio = dummy[line, 'start']
            final = dummy[line, 'end']
            
            # vamos a contar los exones que pertenezcan a ese cromosoma y contengan al bloque
            query = paste('select count (*) from data where seqname="', chr, '" and feature="exon" and start<', toString(inicio), ' and end>', toString(final), sep ='')
            num_exon = dbGetQuery(db2, query)
            
            # guardamos la fila completa para escribirla
            row = dummy[line, ]
            
            # si contamos sólo uno, escribimos en exones; si contamos cero, lo escribimos en intrones; si contamos otra cosa, hay un error. Llevaremos la cuenta de los errores.
            if(num_exon == 1){
                
                exones = rbind(exones, row)
                
            }else if(num_exon == 0){
                
                intrones = rbind(intrones, row)
                
            }else{
                
                num_errores = num_errores + 1
                
            }
            
            
        }
        # recuperamos el valor de expresión para ese gen
        expression = dbGetQuery(db3, paste('select expression from ', tejido, ' where gen="', gen, '"'))
        
        # calculamos los promedios de beta para exones e intrones
        avg_exones = mean(exones[,6])
        avg_intrones = mean(intrones[,6])
        
        # diseñamos las líneas que serán escritas a los archivos temporales
        # Las columnas de las tablas originales son: (closest_gene_id, closest_gene_name,
        # chr, start, end, beta, block_location, closest_gene_strand)
        # Las columnas de la tabla de salida son: (closest_gene_id, closest_gene_name,
        # chr, beta_avg, expression, block_location, closest_gene_strand)
        line_exones = exones[1,1:3]
        line_exones = c(line_exones, avg_exones, expression, 'exon', exones[1,8])
        
        line_intrones = intrones[1,1:3]
        line_intrones = c(line_intrones, avg_intrones, expression, 'intron', intrones[1,8])
        
        # escribimos las líneas a los archivos temporales
        write.table(t(line_exones), file=tmpfile_exones, col.names=FALSE, row.names=FALSE, append=TRUE, sep='\t')
        write.table(t(line_intrones), file=tmpfile_intrones, col.names=FALSE, row.names=FALSE, append=TRUE, sep='\t')
        
        
        print(paste('Gen ', gen, ' de ', tabla, 'terminado.'))
     
    }
    
    # creamos las tablas en la base de datos para exones e intrones
    dbGetQuery(db, paste('create table ', tejido, '_avg_exones(closest_gene_id text, closest_gene_name text, chr text, beta_avg real, expression real, block_location text, closest_gene_strand integer)'))
    dbGetQuery(db, paste('create table ', tejido, '_avg_intrones(closest_gene_id text, closest_gene_name text, chr text, beta_avg real, expression real, block_location text, closest_gene_strand integer)'))
    
    # cargamos los datos desde archivo a la base de datos
    dbWriteTable(db, paste(tejido, '_avg_exones', sep=''), tmpfile_exones, sep='\t', header = FALSE, append = TRUE)
    dbWriteTable(db, paste(tejido, '_avg_intrones', sep=''), tmpfile_intrones, sep='\t', header = FALSE, append = TRUE)
    
    # borramos los archivos temporales
    file.remove(tmpfile_exones)
    file.remove(tmpfile_intrones)
    
    # escribimos en un archivo el número de errores hallados para esa tabla
    errores = file(paste(tabla, '_errores.txt', sep=''))
    write(toString(num_errores), file = errores)
}

#############################################################################################


# checamos la lista de tablas de la base de datos
tablas = dbListTables(db)

# inicializamos la lista de tablas que nos interesan
tablas_originales = c()

# iteramos por todas para ver cuáles son las originales
for(tabla in tablas){
    # guardamos el número de caracteres en el string
    len = nchar(tabla)
    
    # checamos si el nombre de la tabla tiene 8 caracteres o menos
    if(substr(tabla, len-9, len) == 'Distancias'){
        # guardamos el nombre de esa tabla en la lista
        tablas_originales = c(tablas_originales, tabla)
    }
}


# iteramos por la lista de tablas que nos interesa
for(tabla in tablas_originales){
    # le aplicamos la función a esa tabla
    separar_ex_int(tabla)
}
