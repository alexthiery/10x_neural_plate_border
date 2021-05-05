class Helper {
    private static ArrayList globArray(array, glob){
        ArrayList file_array = []
        for(def file:array){
            if(file.toString().matches(glob)){
                file_array.add(file)
            }
        }
        return file_array
    }
}