package com.company;

import java.io.*;
import java.nio.file.*;
import java.nio.file.attribute.BasicFileAttributes;

public class Main {

    public static class PrintFiles extends SimpleFileVisitor<Path> {

        // Print information about
        // each type of file.
        @Override
        public FileVisitResult visitFile(Path file, BasicFileAttributes attr) throws IOException {
            if (attr.isSymbolicLink()) {
                System.out.format("Symbolic link: %s ", file);
            } else if (attr.isRegularFile()) {
                System.out.format("Regular file: %s ", file);

               /*if (file.toString().contains("aligned")) {
                    File f = new File(file.toAbsolutePath().toString());
                    if (f.exists()) {
                        f.delete();
                    }

                    return FileVisitResult.CONTINUE;
                }*/


                String seq1 = "";
                String seq2 = "";

                FileReader fr = null;
                BufferedReader br = null;

                int s1, s2;

                try {
                    fr = new FileReader(file.toAbsolutePath().toString());
                    br = new BufferedReader(fr);

                    String line = "";
                    StringBuilder sb = new StringBuilder();

                    line = br.readLine();

                    while (true) {
                        line = br.readLine();
                        if (line.startsWith(">")) break;
                        sb.append(line);
                    }

                    seq1 = sb.toString();

                    s1 = seq1.length();

                    sb = new StringBuilder();

                    while ((line = br.readLine()) != null) {
                        sb.append(line);
                    }

                    seq2 = sb.toString();

                    s2 = seq2.length();
                } catch(Exception e) {
                    e.printStackTrace();
                } finally {
                    fr.close();
                    br.close();
                }

                FileWriter fw = null;
                BufferedWriter bw = null;

                try {
                    String fileName = file.getFileName().toString().replace(".fasta", "-aligned.fasta");

                    fw = new FileWriter("/home/martin/aligned/" + fileName);
                    bw = new BufferedWriter(fw);

                    bw.write(seq1);
                    bw.write("\n");
                    bw.write(seq2);

                } catch(Exception e) {

                } finally {
                    fw.close();
                    // bw.close();
                }
            } else {
                System.out.format("Other: %s ", file);
            }
            System.out.println("(" + attr.size() + "bytes)");
            return FileVisitResult.CONTINUE;
        }

        // Print each directory visited.
        @Override
        public FileVisitResult postVisitDirectory(Path dir, IOException exc) {
            System.out.format("Directory: %s%n", dir);
            return FileVisitResult.CONTINUE;
        }

        // If there is some error accessing
        // the file, let the user know.
        // If you don't override this method
        // and an error occurs, an IOException
        // is thrown.
        @Override
        public FileVisitResult visitFileFailed(Path file, IOException exc) {
            System.err.println(exc);
            return FileVisitResult.CONTINUE;
        }
    }


    public static void main(String[] args) throws IOException {
        Path startingDir = Paths.get("/home/martin/Source/bioinf_project/database/outputs");
        PrintFiles pf = new PrintFiles();
        Files.walkFileTree(startingDir, pf);
    }
}
