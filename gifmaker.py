def gifmaker(gif_name,figpath):
    import glob
    import os

    file_list = glob.glob(figpath + '*.png') # Get all the pngs in the current directory
    list.sort(file_list, key=lambda x: int(x.split('_')[1].split('.png')[0])) # Sort the images by #, this may need to be tweaked for your use case

    with open('image_list.txt', 'w') as file:
        for item in file_list:
            file.write("%s\n" % item)

    os.system('convert @image_list.txt {}.gif'.format(gif_name)) # On windows convert is 'magick'
    

