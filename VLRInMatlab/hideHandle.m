function hideHandle( handle )
if ishandle( handle )
  set( handle, 'HandleVisibility', 'off' );
end