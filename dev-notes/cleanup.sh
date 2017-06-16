while true; do
    read -p "Remove/Replace all site-related files?" yn
    case $yn in
        [Yy]* )
        rm *.pdf *.html *.css *.json
        rm -rf libs
        mv dev-notes-output/* .;
        break;;
        [Nn]* ) exit;;
        * ) echo "Please answer yes or no.";;
    esac
done
