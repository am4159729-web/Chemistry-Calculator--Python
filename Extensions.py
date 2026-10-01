#I updated Option 7 so it can automatically talk to our new, separate add-ons feature

elif choice == '7':
  if not HISTORY:
    print("\nNo calculations recorded in this session yet.\n")
  else:
    print(f"\n=== SESSION HISTORY ({len(HISTORY)} entries) ===")
    for idx, item in enumerate(HISTORY, start=1):
    print(f"[{idx}] {item.operation} -> Result: {item.result}")
    print("==========================================\n")
                    
   # ADD THESE THREE LINES BELOW YOUR HISTORY PRINT LOOP 
    
    export_choice = input("Would you like to export this history to a Markdown report? (y/n): ").strip().lower()
    if export_choice == 'y':
      import extensions
      extensions.export_history_to_markdown()
